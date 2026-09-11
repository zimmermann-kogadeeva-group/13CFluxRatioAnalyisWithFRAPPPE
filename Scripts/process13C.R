library(tidyverse)
library(readxl)
library(here)
# Functions to process labeling data

# utils

spec <- function(df) {
    df %>% 
        pull(identifier) %>% 
        unique() %>% 
        sort() %>% 
        return()
}

# turn pvalue in proper scientific writing notion

format_pvalue <- function(pvalue) {
  pval_str <- formatC(pvalue, format = "e", digits = 2)
  pval_out <- sub("e", " %*% 10^", pval_str)
  return(pval_out)
  
}

pivot_wider_fortables <- function(df, mode = NA) {
  
    output <- df %>% 
      mutate(sample_name = paste0(experiment, "_", strain, "_", medium, "_", brep)) %>% 
      pivot_wider(id_cols = c("identifier", "isotopologue", "mzCloud", "mzVault", "rt_min"), names_from = "sample_name", values_from = "relab") 
    
    if(!(is.na(mode))) {
        
        output <- output %>% 
            mutate(mode = mode)
      
    }
    
    
    return(output)
}


# 1. general function to read in data from CD3.4 and turn it into a relative abundance table

format_CD34_output <- function(
    path_to_CD34_output,
    path_to_metadata, 
    path_to_annotation_corrections = NA
) { 
    
    # split levels of the dataframe
    
    data <- read_xlsx(here(path_to_CD34_output))
    metadata <- read_csv(path_to_metadata)
    colnames(data)[11] <- "delta_ppm"
    colnames(data)[12] <- "calc_mw"
    colnames(data)[13] <- "rt_min"
    colnames(data)[14] <- "area_max"

    # Fix names of features where the annotation is not 100% correct
    
    if(!(is.na(path_to_annotation_corrections))) {
        
        annotation_corrections <- read_csv(here(path_to_annotation_corrections)) %>% mutate(rt_min = as.character(rt_min))
        
        data <- data %>% 
            left_join(annotation_corrections, by = c("Name", "rt_min")) %>% 
            mutate(Name = ifelse(is.na(new_name), Name, new_name))
         
    }
    
    data2 <- data %>% 
        mutate(level = ifelse(grepl("Polarity|\\+|\\-", `Annot. Source: mzVault Search`), 2, 1)) %>% 
        mutate(identifier = ifelse(level == 1, paste0(Name, "_", calc_mw, "_", rt_min), NA)) %>% # make new id col where name is in when lvl1, then fill name down
        fill(identifier, .direction = "down")
        
    lvl1 <- data2 %>% # names, annotations ect
        filter(level == 1) %>% 
        select( # select relevant features
            identifier, 
            Name, 
            Formula, 
            `Annot. Source: mzCloud Search`, 
            `Annot. Source: mzVault Search`, 
            delta_ppm,
            rt_min, 
            area_max, 
            MS2, 
            Background
        )

    lvl2 <- filter(data2, level == 2)
    colnames(lvl2)[1:(ncol(lvl2) - 2)] <- lvl2[1, 1:(ncol(lvl2) - 2)]        
    
    # create joint dataframe
    df <- lvl2 %>% 
        select( # remove all junk
            contains("Exchange Rate"), 
            identifier,
            "Study File ID", 
            "Area"
        ) %>% 
        pivot_longer(
            cols = contains("Exchange Rate"), 
            names_to = "isotopologue",
            values_to = "relab"
        ) %>% 
        filter(isotopologue != relab) %>% 
        rename("study_file_id" = "Study File ID") %>% 
        left_join(metadata, by = "study_file_id") %>% # bind with other infos
        left_join(lvl1, by = "identifier") 
        
    # rename some columns for ease of useage
    
    df_out <- df %>% 
        mutate(
            isotopologue = gsub(".*: ", "m+", isotopologue), 
            relab = as.numeric(relab) / 100, 
            delta_ppm = as.numeric(delta_ppm), 
            rt_min = as.numeric(rt_min), 
            Area = as.numeric(Area)
        ) %>%
        rename(
            "mzVault" = "Annot. Source: mzVault Search", 
            "mzCloud" = "Annot. Source: mzCloud Search"
        )
    
    return(df_out)
}

# bind two dfs from positive and negative mode

bind_neg_pos <- function(df_neg, df_pos) {
    
    output <- bind_rows(
        df_neg %>% mutate(mode = "neg"), 
        df_pos %>% mutate(mode = "pos")
    )
    
    return(output)
}

# calculate mean relative abundances 

calc_mean_relabs <- function(relabs) {
      
    if(is_grouped_df(relabs) == FALSE) {
        stop("Needs to supply a grouped dataframe, otherwise it does not make sense to apply this function! Exiting...")
    }
    
    mean_relabs <- relabs %>% 
        summarise(
            mean = mean(relab, na.rm = TRUE), 
            sd = sd(relab, na.rm = TRUE)
        )
     
    return(mean_relabs)
     
}

# remove features that when rt is known but outside of tolerance window of 1min and add level information to the features

filter_rtlib <- function(relabs, rts) { 
    
    rts <- rts %>% 
        mutate(in_rtlist = TRUE) %>% 
        mutate(rt_ul = rt_min + 1, rt_ll = rt_min - 1) %>% 
        select(-rt_min)
    
    relabs_noinfo <- relabs %>% 
        filter(!(Name %in% rts$compound_name)) %>%  # filter out those features that are not part of the library, so on those we cannot make any claims
        mutate(in_rtlist = FALSE)
        
    relabs_rtsinfo <- relabs %>% 
        filter(Name %in% rts$compound_name) %>% 
        left_join(rts, by = join_by(Name == compound_name, between(rt_min, rt_ll, rt_ul))) %>% 
        filter(in_rtlist)
    
    # join together
    
    relabs_filt <- bind_rows(
        relabs_noinfo, 
        relabs_rtsinfo
    )
    
    return(relabs_filt)
}

# function to figure out which compounds are level 1 or level 2

add_mslvl <- function(relabs) {
    
    if(!("in_rtlist" %in% colnames(relabs))) {
        stop("no column 'in_rtlist' found! Please first run the function 'filter_rtlib' to add information about retention time accuracy to the data! Exiting...")
    }
    
    relabs_wlvl <- relabs %>% 
        mutate(mslvl = ifelse(mzVault == "Full match" & in_rtlist, 1, NA)) %>% 
        mutate(mslvl = ifelse(is.na(mslvl) & (mzVault == "Full match" | mzCloud == "Full match"), 2, mslvl))
        
    return(relabs_wlvl)
        
}

# remove features with no MS2 annotation

filter_MS2 <- function(relabs) {
    
    relabs_filt <- relabs %>% 
        filter(mzVault == "Full match" | mzCloud == "Full match")
    
    return(relabs_filt)
}

# remove features which do not surpass the background level

filter_background <- function(relabs, background_multiplier) {
    
    blanks <- relabs %>%
        filter(strain == "VBDBT") %>% 
        select(identifier, experiment, Area, study_file_id, brep) %>% 
        unique() %>% 
        group_by(identifier, experiment) %>% 
        summarise(background_mean = mean(Area))

    relabs_brem <- relabs %>% 
        left_join(blanks, by = c("identifier", "experiment")) %>% 
        mutate(background_mean = ifelse(is.na(background_mean), 10000, background_mean)) %>% 
        mutate(is_above_background = ifelse(Area >= background_multiplier*background_mean, TRUE, FALSE)) %>% 
        filter(is_above_background)
        
    return(relabs_brem)

}

# remove features not observed in a specific number of replicates. provide 

filter_nobs <- function(relabs, nobs) {
        
    if(is_grouped_df(relabs) == FALSE) {
        stop("Needs to supply a grouped dataframe, otherwise it does not make sense to apply this function! Exiting...")
    }
    
    relabs_nobs <- relabs %>% 
        mutate(identified_in = length(unique(brep))) %>% 
        filter(identified_in == nobs) %>% 
    
    return(relabs_nobs)
        
}

filter_sd <- function(relabs, sd_cutoff, group_by_medium = FALSE, group_by_strain = FALSE) {
    
    mean_relabs <- calc_mean_relabs(relabs %>% group_by(strain, medium, experiment, isotopologue, identifier))
        
    mean_relabs_to_remove <- mean_relabs %>% 
        ungroup()
        
    # decide on the grouping criterion for the  standard deviation criterion
    
    if(group_by_medium & !group_by_strain) {
        print("group by medium")
        mean_relabs_to_remove <- group_by(mean_relabs_to_remove, identifier, medium)
    } else if (group_by_strain & !group_by_medium) {
        print("group by strain")
       mean_relabs_to_remove <- group_by(mean_relabs_to_remove, identifier, strain)
    } else if (group_by_strain & group_by_medium) {
        print("group by medium and strain")
        mean_relabs_to_remove <- group_by(mean_relabs_to_remove, identifier, strain, medium)
    } else {
        print("only identifier grouping")
        mean_relabs_to_remove <- group_by(mean_relabs_to_remove, identifier)
    }
    
    mean_relabs_to_remove <- mean_relabs_to_remove %>% 
        mutate(max_sd = max(sd, na.rm = TRUE)) %>% 
        filter(max_sd >= sd_cutoff) %>% 
        select(medium, identifier) %>% 
        mutate(toremove = TRUE) %>% 
        unique()
    
    relabs_filt <- relabs %>% 
        left_join(mean_relabs_to_remove) %>% 
        filter(is.na(toremove))

    return(relabs_filt)            
    
}

# mark features occuring multiple times in the same method

filter_not_in_all <- function(relabs) {
  
  mean_relabs <- calc_mean_relabs(relabs %>% group_by(strain, medium, experiment, isotopologue, identifier))
  
  groups_medium <- mean_relabs %>% pull(medium) %>% unique()
  groups_strain <- mean_relabs %>% pull(strain) %>% unique()
  
  mean_relabs <- mean_relabs %>% 
    group_by(identifier, isotopologue) %>% 
    filter(n() == (length(groups_medium) * length(groups_strain)))
  
  relabs_filt <- filter(relabs, identifier %in% mean_relabs$identifier)
  
  return(relabs_filt)
  
}

filter_m0_threshold <- function(relabs, m0_group, m0_threshold) {
    
    features_qc_filtered <- relabs %>%
        filter(medium == m0_group & isotopologue == "m+0" & relab >= m0_threshold) %>%
        pull(identifier) %>%
        unique()
    
    relabs_filt <- filter(relabs, !(identifier %in% features_qc_filtered))
    
    return(relabs_filt)
    
}

identify_duplicate_feats <- function(relabs) {
  
  duplicate_feats <- relabs %>% 
    select(Name, rt_min, identifier, area_max) %>% 
    unique() %>% 
    group_by(Name) %>% 
    mutate(marked_duplicate = ifelse(max(area_max) == area_max, FALSE, TRUE)) %>% 
    ungroup() %>% 
    select(identifier, marked_duplicate)
  
  relabs <- relabs %>% 
    left_join(duplicate_feats, by = "identifier") %>% 
    return()
  
}

# calculate relative amount of labeled carbon

calc_perc_lab_carb <- function(relabs) {

    mean_relabs <- calc_mean_relabs(relabs %>% group_by(strain, medium, experiment, isotopologue, identifier))    

    perc_lab_carb <- mean_relabs %>% 
        mutate(ncarb = as.numeric(gsub("m+", "", isotopologue))) %>% 
        group_by(strain, identifier, medium) %>% 
        summarise(
            perc_labeled = (sum(ncarb*mean)/(n()-1))
        )
    
    return(perc_lab_carb)
}

calc_cosine_similarity <- function(relabs, ref_strain = NA, ref_medium = NA) {
    
    cosim <- function(A,B) { (sum(A*B))/sqrt((sum(A^2))*(sum(B^2))) }
    
    mean_relabs <- calc_mean_relabs(relabs %>% group_by(strain, medium, experiment, isotopologue, identifier))
    
    if(!(is.na(ref_strain)) & is.na(ref_medium)) {
        
        ref_group_df <- mean_relabs %>%     
            filter(strain == ref_strain) %>% 
            rename("ref_mean" = "mean") %>% 
            ungroup() %>%
            select(medium, identifier, isotopologue, ref_mean)
        
        cosim <- mean_relabs %>% 
            left_join(ref_group_df, by = c("medium", "isotopologue", "identifier")) %>% 
            group_by(strain, medium, identifier) %>% 
            summarise(
                cosim = cosim(mean, ref_mean)
            )
    
    } else if(!(is.na(ref_medium)) & is.na(ref_strain)) { 
        
        ref_group_df <- mean_relabs %>% 
            filter(medium == ref_medium)%>% 
            rename("ref_mean" = "mean") %>% 
            ungroup() %>% 
            select(strain, identifier, isotopologue, ref_mean)
        
        cosim <- mean_relabs %>% 
            left_join(ref_group_df, by = c("strain", "isotopologue", "identifier")) %>% 
            group_by(strain, medium, identifier) %>% 
            summarise(
                cosim = cosim(mean, ref_mean)
            )
    
    } else if(!(is.na(ref_medium)) & !(is.na(ref_strain))) { 
    
        ref_group_df <- mean_relabs %>% 
            filter(medium == ref_medium & strain == ref_strain) %>% 
            rename("ref_mean" = "mean") %>% 
            ungroup() %>% 
            select(identifier, isotopologue, ref_mean)
    
        cosim <- mean_relabs %>% 
            left_join(ref_group_df, by = c("isotopologue", "identifier")) %>% 
            group_by(strain, medium, identifier) %>% 
            summarise(
                cosim = cosim(mean, ref_mean)
            )
    
    } else {
        stop("Pick a reference strain or medium for comparison! Exiting...")
    }
}

# function to extract labeling patterns for modeling

extract_patterns_for_modeling <- function(relabs, path_to_index) {
    
    index <- read_csv(path_to_index) %>% 
      select(-mode)
    
    if(!("mode" %in% colnames(relabs))) {
        stop("First join positive and negative dataframes together! Exiting...")
    }
    
    mean_relabs <- calc_mean_relabs(relabs %>% group_by(strain, medium, experiment, isotopologue, identifier, mode)) %>% # bind names and index for model mapping
      left_join(relabs %>% select(identifier, Name) %>% unique(), by = "identifier") %>% 
      inner_join(index, by = join_by("Name" == "metabolite_data"))
    
    # adapt names to match names in modeling framework
    
    output <- mean_relabs %>%
        ungroup()  %>% 
        mutate(ncarb = as.numeric(gsub("m+", "", isotopologue))) %>% 
        group_by(identifier) %>% 
        mutate(nmax = max(ncarb)) %>% 
        mutate(substr = map_chr(nmax, ~ paste0(1:.x, collapse = ""))) %>% 
        mutate(isot_model_name = paste0(metabolite_model, "_", substr, "_", ncarb)) %>%        # redesign name
        filter(!(is.na(metabolite_model))) %>% 
        ungroup() %>% 
        select(identifier, strain, medium, isot_model_name, mean, sd) 
    
    return(output)
        
}


# plotting functions

plot_labeling_pattern <- function(relabs, ident, output_dir, by_strain = FALSE, by_medium = FALSE) {
    
    if(!by_strain & !by_medium) {
        stop("Choose a plot grouping, either by strain or by medium! Exiting...")
    }
    
    relabs <- mutate(relabs, isotopologue = factor(isotopologue, levels = paste0("m+", c(0:40)), labels = paste0("m+", c(0:40))))  
    
    relabs_ident <- filter(relabs, identifier == ident)
    
    mean_relabs <- calc_mean_relabs(relabs_ident %>% group_by(strain, medium, experiment, isotopologue, identifier))
    
    f <- gsub("\\/", "_", ident) # prevent name containing backslash from messing up the plot path
    
    lvl <- as.character(unique(pull(relabs_ident, mslvl)))
    
    if(by_strain) {
        
        colors <- c(
            "Buni" = "#ef476f", 
            "Ecoli" = "#118ab2",
            "Pvul" = "#06d6a0"         
        )
        
        p <- ggplot(data = mean_relabs, mapping = aes(x = isotopologue, y = mean, fill = strain))
        
    } else {
            
        colors <- c(
            Glucose   = "#56B4E9",
            Ribose    = "#CC79A7",
            Thymidine = "#D55E00",
            Adenosine = "#009E73",
            Uridine   = "#E69F00",
            Cytidine  = "#0072B2"
        )
        
        p <- ggplot(data = mean_relabs, mapping = aes(x = isotopologue, y = mean, fill = medium))
    
    }

    p <- p + 
        geom_col(position="dodge") + 
        geom_errorbar(mapping = aes(x = isotopologue, y = mean, ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9), color = "black") + 
        geom_point(data = relabs_ident, mapping = aes(x = isotopologue, y = relab), position = position_dodge(width = 0.9)) + 
        xlab("") + 
        ylab("mean relative \nisotopologue abundance") + 
        theme_bw() + 
        scale_fill_manual(values = colors) + 
        ylim(0,1) + 
        facet_wrap(. ~ identifier) + 
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            strip.background = element_rect(fill = "white")
        )
    

    ggsave(
        p, 
        filename = paste0(f, "_", lvl, ".pdf"), 
        path = output_dir,
        width = 4.5,
        height = 2.5
    )   
}

plot_labeling_pattern_timeseries <- function(relabs, ident, output_dir) {
    
    relabs <- mutate(relabs, isotopologue = factor(isotopologue, levels = paste0("m+", c(0:40)), labels = paste0("m+", c(0:40)))) %>% 
        mutate(timepoint = factor(timepoint, levels = c("0", "10", '20', "40", "80", "160", "320")))
    
    relabs_ident <- filter(relabs, identifier == ident)
    
    mean_relabs <- calc_mean_relabs(relabs_ident %>% group_by(strain, medium, experiment, isotopologue, identifier, timepoint))
    
    f <- gsub("\\/", "_", ident) # prevent name containing backslash from messing up the plot path
    
    lvl <- as.character(unique(pull(relabs_ident, mslvl)))

    colors <- c(
        "m+0" = "#e4e7e4", 
        "m+1" = "#ff595e", 
        "m+2" = "#8ac926",
        "m+3" = "#ffca3a",
        "m+4" = "#1982c4",
        "m+5" = "#6a4c93",
        "m+6" = "#0a1647"
    )    
    
     
    p <- ggplot(data = mean_relabs, mapping = aes(x = timepoint, y = mean, fill = isotopologue)) + 
        geom_bar(position="stack", stat = "identity")  + 
        xlab("") + 
        ylab("mean relative \nisotopologue abundance") + 
        theme_bw() + 
        scale_fill_manual(values = colors) + 
        # ylim(0,1) + 
        facet_wrap(. ~ identifier) + 
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            strip.background = element_rect(fill = "white"), 
            strip.text.x = element_text(size = 8)
        ) + 
        xlab("time [s]")

    ggsave(
        p,
        filename = paste0(f, ".pdf"), 
        path = output_dir,
        width = 4.5,
        height = 2.5
    )
}

plot_perc_carb <- function(df, ref = NA, comp = NA, title, output_dir) {
    
    if(ref %in% c("Buni", "Pvul", "Ecoli")) { # ensure proper axis labels of species strings
      
      ital = TRUE
      
      if (ref == "Buni") {
        ref_str = "B. uniformis"
      } else if(ref == "Pvul") {
        ref_str = "P. vulgatus"
      } else if(ref == "Ecoli") {
        ref_str = "E. coli"
      }
      
      if (comp == "Buni") {
        comp_str = "B. uniformis"
      } else if(comp == "Pvul") {
        comp_str = "P. vulgatus"
      } else if(comp == "Ecoli") {
        comp_str = "E. coli"
      }
      
    } else {
      ref_str <- ref
      comp_str <- comp
      ital = FALSE
    }
      
    refs <- df %>% 
        filter(strain == ref | medium == ref) %>% 
        ungroup() %>% 
        select(identifier, perc_labeled) %>% 
        rename("perc_ref" = "perc_labeled")
    
    df2 <- df %>% 
        filter(strain == comp | medium == comp) %>% 
        left_join(refs, by = c("identifier")) %>% 
        filter(!(is.na(perc_ref)))
    
    # calculate number of features and significance of depletion
    
    nfeat <- nrow(df2)
    
    pval <- wilcox.test(
      df2$perc_labeled,
      df2$perc_ref, 
      alternative = "two.sided", 
      paired = TRUE
    )$p.value
    
    pvalue_formatted <- format_pvalue(pval)
    
    
    if(ital) { # turn strings italic in case it is a species name
      
      xtitle <- bquote(atop("fraction of labeled carbon in", italic(.(ref_str)) ~ " metabolites"))
      ytitle <- bquote(atop("fraction of labeled carbon in", italic(.(comp_str)) ~ " metabolites"))
      
    } else {
      
      xtitle <- bquote(atop("fraction of labeled carbon in", .(ref_str) ~ " condition metabolites"))
      ytitle <- bquote(atop("fraction of labeled carbon in", .(comp_str) ~ " condition metabolites"))
      
    }
    
    ggplot(df2, aes(x = perc_ref, y = perc_labeled)) + 
        geom_point(alpha = 0.65) + 
        annotate("text", x = 0, y = 0.95,
                 label = paste0("'n = ", nfeat, ", ' ~ italic(p) ~ ' = ' ~ ", pvalue_formatted),
                 size = 4, parse = TRUE, hjust = 0) +
        labs(x = xtitle, y = ytitle) + 
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
        theme_bw() + 
        facet_wrap(medium ~ .) + 
        ylim(0,1) + 
        xlim(0,1) + 
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.text.x = element_text(size = 11), 
            axis.text.y = element_text(size = 11), 
            strip.text = element_text(size = 13), 
            legend.text = element_text(size = 11), 
            legend.title = element_text(size = 13), 
            axis.title.x = element_text(size = 13), 
            axis.title.y = element_text(size = 13), 
            strip.background = element_rect(fill = "white")
        )  
    
    ggsave(
        filename = title, 
        path = output_dir,
        width = 3.5,
        height = 3.5
    )
}

# plot boxplot 

plot_cosine_similarities_boxplot <- function(cosine_sim, comp, title, output_dir) { 
    
    cosine_sim_filt <- cosine_sim %>% 
      filter(!(is.na(cosim))) 
    
    nfeat <- cosine_sim_filt %>% 
      group_by(!! sym(comp)) %>% 
      tally()
    
    ggplot(cosine_sim_filt, aes(x = !! sym(comp), y = cosim)) + 
        geom_boxplot(outlier.shape = NA) + 
        geom_point(position = position_jitter(width = 0.2), alpha = 0.6) + 
        geom_text(data = nfeat, aes(x = !! sym(comp), y = 1, label = paste0("n=", n), vjust = -1), 
                  inherit.aes = FALSE, size = 4) +
        theme_bw() + 
        facet_wrap(strain ~ .) + 
        coord_cartesian(ylim = c(0, 1.1), clip = "off") + 
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            strip.background = element_rect(fill = "white"), 
            axis.text.x = element_text(size = 11, angle = 30, hjust = 1), 
            axis.text.y = element_text(size = 11), 
            strip.text = element_text(size = 13), 
            legend.text = element_text(size = 11), 
            legend.title = element_text(size = 13), 
            axis.title.x = element_text(size = 13), 
            axis.title.y = element_text(size = 13),
            plot.margin = margin(t = 20, r = 5, b = 5, l = 5)
        ) + 
        xlab("") + 
        ylab("Cosine similarity")

        
    ggsave(
        filename = title, 
        path = output_dir,
        width = 4.5,
        height = 3.5
    )
}

plot_cosine_similarities_scatterplot <- function(df, x, y, title, output_dir, x_str, y_str) { 
    
    if(x %in% c("Buni", "Pvul")) { # ensure proper axis labels of species strings
      
        if (x == "Buni") {
          x_str = "B. uniformis"
        } else if(y == "Pvul") {
          x_str = "P. vulgatus"
        } 
        
        if (y == "Buni") {
          y_str = "B. uniformis"
        } else if(y == "Pvul") {
          y_str = "P. vulgatus"
        } 
        
    }
  
    cosine_sim_x <- df %>% 
        filter(strain == x | medium == x) %>% 
        ungroup() %>% 
        select(identifier, cosim, medium) %>% 
        rename("cosim_x" = "cosim")
        
    cosine_sim_y <- df %>% 
        filter(strain == y | medium == y) %>% 
        ungroup() %>% 
        select(identifier, cosim, medium) %>% 
        rename("cosim_y" = "cosim")
        
    cosine_joint <- full_join(cosine_sim_x, cosine_sim_y, by = c("identifier", "medium")) %>% filter(!(is.na(cosim_x))) %>% filter(!(is.na(cosim_y)))
    nfeat <- as.character(nrow(cosine_joint))
    
    xtitle <- bquote(atop("Cosine similarity between", italic(.(x_str)) ~ " and " ~ italic("E. coli") ~ " MDVs"))
    ytitle <- bquote(atop("Cosine similarity between", italic(.(y_str)) ~ " and " ~ italic("E. coli") ~ " MDVs"))
    
    ggplot(cosine_joint, aes(cosim_x, cosim_y)) + 
        geom_point(alpha = 0.75) + 
        annotate("text", x = 0, y = 0.95,
                 label = paste0("n = ", nfeat),
                 size = 5, hjust = 0) + 
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
        labs(x = xtitle, y = ytitle) + 
        facet_wrap(medium ~ .) + 
        theme_bw() + 
        ylim(0,1) + 
        xlim(0,1) + 
        scale_color_manual(values = colors) + 
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 11), 
            axis.text.y = element_text(size = 11), 
            strip.text = element_text(size = 13), 
            legend.text = element_text(size = 11), 
            legend.title = element_text(size = 13), 
            axis.title.x = element_text(size = 13), 
            axis.title.y = element_text(size = 13),
            strip.background = element_rect(fill = "white")
        )

    ggsave(
        filename = title, 
        path = output_dir,
        width = 3.5,
        height = 3.5
    )
        
}


