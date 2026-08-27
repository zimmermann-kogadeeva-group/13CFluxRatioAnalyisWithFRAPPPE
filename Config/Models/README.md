
# Models

## Toy

- `Toy/freeflux.tsv`
    - model definition taken from one of the examples in freeflux package repo
    - Used to test how freeflux works and how to use these model defitions with
    hopsy to generate fluxes

## Ecoli

- `Ecoli/reactions_model_inca_ecoli.csv`
    - model to compare hopsy/freeflux with sumoflux 
    - from `/g/mazimmer/mazimmer/Projects/MZK009_Past/MZK009A_SUMOFLUX/MZK009Aa_SUMOFLUX_cofeed/Data/Ecoli_data_for_SUMOFLUXpy_testing/reactions_model_inca_ecoli.csv`

- `Ecoli/reactions_model_inca_ecoli.tsv`
    - Based on `Ecoli/reactions_model_inca_ecoli.csv`, but tidied up - had to remove extra commas

- `Ecoli/reactions_model_inca_ecoli_corrected.tsv`
    - Based on `Ecoli/reactions_model_inca_ecoli.tsv`, but fixed reversibility of some reactions (fba, bpg, eno, gapdh)
    - Used in `Notebooks/troubleshooting_ecoli.ipynb`

- `Ecoli/freeflux_syn.xlsx`
    - Model from freeflux github repo
    - Used in `Notebooks/ecoli_example.ipynb`

- `Ecoli/freeflux_w_co2_influx.tsv`
    - Based on `Ecoli/freeflux_syn.xslx`
    - Added co2 influx with co2in reaction

- `Ecoli/freeflux_exp.xlsx`
    - Alternative model from freeflux that is in the `experimenal_data`
    subfolder in their repo
    - Used in `Notebooks/ecoli_example.ipynb`
    - For comparison with syn counterpart

- `Ecoli/Ecoli_individualBiomass.tsv`
    - Model for species comarison of Ecoli Buni and Pvul
    - Used in DBT005G17_evaluation_Ecoli_Buni_Pvul.ipynb
    
    
## Buniformis

- `Buniformis/Buniformis_biomass.tsv`
    - baseline model for B. uniformis
    - Derived from the baseline model formulated in DBT005G1 (scratch_3) and further refined in DBT005G4
    - Final refinement conducted in DBT005G10 (interim models not uploaded for now)
    
- `Buniformis/Buniformis_individualBiomass.tsv`
    - fixed baseline model for B. uniformis, dreived from Buniformis_biomass.tsv
    - has the joint biomass reaction separated into single biomass reaction
    - removed reaction gly2 which turned Threonine into Glycine and Acetyl-CoA, which caused cycling through AA synthesis pathways and m+0 enrichment
    
- `Buniformis/reactions_Buni_marias_import.xlsx`
    - Model translation from the TAC3 model

- `Buniformis/reactions_joint_network.xlsx`
    - Model as reactions_Buni_marias_import.xlsx
    - renamed Ecoli amphibolic reactions by addind a ec tag at the end
    - added amphibolic Buni reactions
    
    
