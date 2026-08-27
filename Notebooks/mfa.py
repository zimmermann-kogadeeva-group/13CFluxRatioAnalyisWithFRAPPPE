import multiprocessing as mp
import os
import re
from functools import wraps
from itertools import permutations, product

import hopsy
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import KFold, cross_val_score, train_test_split
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm


from freeflux import Metabolite, Model, Reaction


def conditional_memoize(env_var_name):
    """
    A decorator factory that conditionally applies diskcache.memoize
    based on the presence of an environment variable.
    """

    def decorator(func):
        try:
            from diskcache import Cache
        except ImportError:
            return func
        else:
            if env_var_name in os.environ:
                # Apply diskcache.memoize if the env var exists
                cache = Cache(os.getenv(env_var_name))

                @wraps(func)
                def wrapped(*args, **kwargs):
                    return cache.memoize()(func)(*args, **kwargs)

                return wrapped

        # Return the original function if the env var does not exist
        return func

    return decorator


def get_model(filename, model_name="unnamed"):
    model = Model(model_name)
    model.read_from_file(filename)
    return model


def make_bounds(stoich_matrix, lb_val, ub_val, specific=None):
    assert lb_val <= ub_val
    lb = lb_val * np.ones(stoich_matrix.shape[1])
    ub = ub_val * np.ones(stoich_matrix.shape[1])
    if specific is not None:
        for r_id, (lb_val, ub_val) in specific.items():
            r_idx = stoich_matrix.columns.get_loc(r_id)
            lb[r_idx] = lb_val
            ub[r_idx] = ub_val

    message = "All upper bounds should be greater than the lower bounds"
    assert np.all(lb <= ub), message

    return lb, ub


def get_xchg_fluxes(fluxes):
    rxs = fluxes.index

    reversible_mask = rxs.str.endswith(("_f", "_b"))
    reversible_rxs = rxs[reversible_mask].str.removesuffix(("_f", "_b")).unique()

    return pd.DataFrame.from_dict(
        {
            r: np.power(
                fluxes.loc[r + "_b"] / fluxes.loc[r + "_f"],
                np.sign(fluxes.loc[r + "_f"] - fluxes.loc[r + "_b"]),
            )
            for r in reversible_rxs
        },
        orient="index",
    )


def get_net_fluxes(fluxes, keep_irreversible=True):
    rxs = fluxes.index

    reversible_mask = rxs.str.endswith(("_f", "_b"))
    reversible_rxs = rxs[reversible_mask].str.removesuffix(("_f", "_b")).unique()
    irreversible_rxs = rxs[~reversible_mask]

    net_fluxes = pd.DataFrame.from_dict(
        {r: fluxes.loc[r + "_f"] - fluxes.loc[r + "_b"] for r in reversible_rxs},
        orient="index",
    )

    if keep_irreversible:
        return pd.concat([net_fluxes, fluxes.loc[irreversible_rxs]])
    else:
        return net_fluxes


def get_reconex_fluxes(fluxes, factor=1, seed=42):
    rxs = fluxes.index

    reversible_mask = rxs.str.endswith(("_f", "_b"))
    reversible_rxs = rxs[reversible_mask].str.removesuffix(("_f", "_b")).unique()
    irreversible_rxs = rxs[~reversible_mask]

    net_fluxes = pd.DataFrame.from_dict(
        {r: fluxes.loc[r + "_f"] - fluxes.loc[r + "_b"] for r in reversible_rxs},
        orient="index",
    )

    rng = np.random.default_rng(seed)
    new_xchg = net_fluxes.apply(
        lambda x: factor * np.sign(x) * rng.uniform(0, np.abs(x)), axis=1
    )

    return pd.concat(
        [
            (net_fluxes + new_xchg)
            .where(net_fluxes > 0, -new_xchg)
            .rename(index=lambda x: x + "_f"),
            new_xchg.where(net_fluxes > 0, -(net_fluxes + new_xchg)).rename(
                index=lambda x: x + "_b"
            ),
            fluxes.loc[irreversible_rxs],
        ]
    ).reindex(rxs)


def reac2array(stoich_mat, reaction):
    return (stoich_mat.columns == reaction).astype(int)


def create_add_constraints_single(stoich_mat, fr_nom, fr_denom, a=0, b=1):
    
    set1 = set(fr_nom)
    set2 = set(fr_denom) - set1

    return np.array([
        np.sum(
            [(1 - b) * stoich_mat.pipe(reac2array, r) for r in set1] + 
            [-b * stoich_mat.pipe(reac2array, r) for r in set2], 
            axis=0
        ),
        np.sum(
            [(a - 1) * stoich_mat.pipe(reac2array, r) for r in set1] + 
            [a * stoich_mat.pipe(reac2array, r) for r in set2], 
            axis=0
        )
    ])


def create_add_constraints(stoich_mat, bounds):
    mat_a = np.concatenate([
        create_add_constraints_single(
            stoich_mat, 
            *ratio, 
            *ratio_bounds
        )
        for ratio_name, (ratio, ratio_bounds) in bounds.items()
    ], axis=0)
    return mat_a


# TODO: replace lb, ub arguments with bounds which is a dictonary of
# arguments for make_bounds func, such that sample_fluxes calls make_bounds inside it
def sample_fluxes(
    stoich_matrix,
    n_samples,
    lb=0,
    ub=10,
    *,
    thinning=10,
    seed=42,
    ss_exclude=None,
    add_constraints=None,
    add_constraints_b=None,
    reconex=None,
    **kwargs,
):
    # Initialize a hopsy problem with only upper bounds and lower bounds, and
    # add the equality constraints afterwards (avoids a bug acc to bartosz)
    if isinstance(lb, (int, float)):
        lb = lb * np.ones(stoich_matrix.shape[1])
    if isinstance(ub, (int, float)):
        ub = ub * np.ones(stoich_matrix.shape[1])

    A = np.concatenate(
        [-np.identity(stoich_matrix.shape[1]), np.identity(stoich_matrix.shape[1])]
    )
    b = np.concatenate([-lb, ub])

    if ss_exclude is not None:
        stoich_matrix = stoich_matrix[~stoich_matrix.index.isin(ss_exclude)]

    if add_constraints is not None:
        add_constraints_b = add_constraints_b or np.zeros(add_constraints.shape[0])
        A = np.concatenate([A, add_constraints])
        b = np.concatenate([b, add_constraints_b])

    problem = hopsy.add_equality_constraints(
        hopsy.Problem(A, b),
        A_eq=stoich_matrix.to_numpy(),
        b_eq=np.zeros(stoich_matrix.shape[0]),
    )

    starting_point = hopsy.compute_chebyshev_center(problem)

    chain = hopsy.MarkovChain(problem, starting_point=starting_point)
    rng = hopsy.RandomNumberGenerator(seed=seed)

    accrate, samples = hopsy.sample(
        chain, rng, n_samples=n_samples, thinning=thinning, **kwargs
    )

    samples = pd.DataFrame(
        samples[0].transpose(),
        index=stoich_matrix.columns,
        columns=[f"s{i}" for i in range(1, samples.shape[1] + 1)],
    )
    if reconex is not None:
        samples = samples.pipe(get_reconex_fluxes, factor=reconex, seed=seed)

    return samples


def get_single_sim_result(
    model_file, fluxes, target_emus, labeling_strategy, sim_name="demo"
):
    model = Model(sim_name)
    model.read_from_file(model_file)

    sim = model.simulator("ss")

    sim.set_target_EMUs(target_emus)
    sim.set_labeling_strategy(**labeling_strategy)
    for r_name, flux in fluxes.items():
        sim.set_flux(r_name, flux)

    sim.prepare()
    res = sim.simulate()
    return res


def get_mdv_dict(mdv_res):
    return {
        f"{emu}_{idx}": x
        for emu in mdv_res.simulated_EMUs
        for idx, x in enumerate(mdv_res.simulated_MDV(emu).value)
    }


def get_mdv_results(
    model_file,
    fluxes_df,
    target_emus,
    labeling_strategy,
    parallel=None,
    as_df=True,
    output_file=None,
):
    if output_file is not None and os.path.isfile(output_file):
        df_freeflux_mdvs = pd.read_csv(output_file, index_col=0)
        if as_df:
            return df_freeflux_mdvs
        else:
            return df_freeflux_mdvs.values

    if parallel is not None:
        vals = [
            (model_file, fluxes, target_emus, labeling_strategy, sim_name)
            for sim_name, fluxes in fluxes_df.to_dict().items()
        ]
        pool = mp.Pool(parallel)
        all_mdv_res = pool.starmap(get_single_sim_result, vals)
    else:
        all_mdv_res = [
            get_single_sim_result(
                model_file, fluxes, target_emus, labeling_strategy, sim_name
            )
            for sim_name, fluxes in tqdm(fluxes_df.to_dict().items())
        ]

    df_freeflux_mdvs = pd.DataFrame(
        [get_mdv_dict(x) for x in all_mdv_res],
        index=[f"s{i:03}" for i in range(1, len(all_mdv_res) + 1)],
    ).transpose()
    if output_file is not None:
        df_freeflux_mdvs.to_csv(output_file)

    if as_df:
        return df_freeflux_mdvs
    else:
        return df_freeflux_mdvs.values


# calculate flux ratios given a dictionary of flux ratios and the corresponding formulas
def get_flux_ratio(data, nominator, denominator, name="flux_ratio"):
    nom_fluxes = []
    denom_fluxes = []
    for r_id in nominator:
        if r_id + "_f" in data.index and r_id + "_b" in data.index:
            nom_fluxes.append(data.loc[r_id + "_f"] - data.loc[r_id + "_b"])
        else:
            nom_fluxes.append(data.loc[r_id])
    for r_id in denominator:
        if (
            r_id + "_f" in data.index and r_id + "_b" in data.index
        ):  # .str.startswith(r_id).sum() == 2:
            denom_fluxes.append(data.loc[r_id + "_f"] - data.loc[r_id + "_b"])
        else:
            denom_fluxes.append(data.loc[r_id])

    return pd.DataFrame(
        np.sum(nom_fluxes, axis=0) / np.sum(denom_fluxes, axis=0), columns=[name]
    )


def get_all_flux_ratios(fluxes, flux_ratios_dict):
    return pd.concat(
        [
            fluxes.pipe(get_flux_ratio, *reac_ratio)
            .assign(flux_ratio_name=name)
            .reset_index(names="iteration")
            for name, reac_ratio in flux_ratios_dict.items()
        ],
        ignore_index=True,
    )


def sample_if_enough(group, samples_per_bin, random_state=87):
    if len(group) >= samples_per_bin:
        return group.sample(samples_per_bin, replace=False, random_state=random_state)
    else:
        return group


def sample_flux_ratio_binned(
    fluxes,
    nominator,
    denominator,
    bins=None,
    samples_per_bin=250,
    random_state=87,
):
    bins = bins or np.linspace(0, 1, 11)
    return (
        get_flux_ratio(fluxes, nominator, denominator)
        .assign(
            binned=lambda x: pd.cut(
                x.flux_ratio, bins=bins, include_lowest=True, right=False
            )
        )
        .reset_index(names="iteration")
        .groupby("binned", as_index=False, group_keys=False)
        .apply(lambda group: sample_if_enough(group, samples_per_bin, random_state))
    )


def subset_fluxes_by_binned_flux_ratio(
    fluxes,
    nominator,
    denominator,
    bins=None,
    samples_per_bin=250,
    random_state=87,
):
    idxs = sample_flux_ratio_binned(
        fluxes, nominator, denominator, bins, samples_per_bin, random_state
    ).iteration
    return fluxes.iloc[:, idxs]


def subset_mdvs_by_binned_flux_ratio(
    mdvs,
    fluxes,
    nominator,
    denominator,
    bins=None,
    samples_per_bin=250,
    random_state=87,
):
    idxs = sample_flux_ratio_binned(
        fluxes, nominator, denominator, bins, samples_per_bin, random_state
    ).iteration
    return mdvs.iloc[:, idxs]


def sample_all_flux_ratio_binned(
    fluxes, flux_ratio_dict, bins=None, samples_per_bin=250, random_state=87
):
    return pd.concat(
        [
            sample_flux_ratio_binned(
                fluxes, *frac, bins, samples_per_bin, random_state
            ).assign(flux_ratio_name=ratio_name)
            for ratio_name, frac in flux_ratio_dict.items()
        ],
        ignore_index=True,
    )


def train_rf(fluxes, mdvs, fr_name, fr_nom, fr_denom, test_size=0.3, random_state=42):
    mdv_train, mdv_test, ratios_train, ratios_test = train_test_split(
        mdvs.transpose(),
        fluxes.pipe(get_flux_ratio, fr_nom, fr_denom).flux_ratio,
        test_size=test_size,
        random_state=random_state,
    )

    kfolds = KFold(n_splits=5, shuffle=True, random_state=random_state)
    rf = RandomForestRegressor(random_state=random_state)

    scores = cross_val_score(
        rf,
        mdv_train,
        ratios_train,
        cv=kfolds,
        scoring="neg_mean_absolute_error",
    )

    score_mean = np.mean(scores)
    score_sd = np.std(scores)

    rf.fit(mdv_train, ratios_train)

    ratios_pred = rf.predict(mdv_test)
    mae = mean_absolute_error(ratios_pred, ratios_test)
    scores = {"mae": mae, "cv_mae_mean": score_mean, "cv_mae_sd": score_sd}
    df_pred_sim = pd.DataFrame({"prediction": ratios_pred, "simulation": ratios_test})
    return rf, scores, df_pred_sim


# Function to replace the standard glucose uptake reaction with a custom uptake reaction
# The metabolite needs to be part of the network already
# By definition, since it is just an uptake reaction, we can use the same carbs
# definition string for both extracellular and intracellular
def replace_uptake_reaction(model, s: str, p: str, carbs: str):

    Subs = Metabolite(s, atoms=carbs)
    Prod = Metabolite(p, atoms=carbs)
    Reac = Reaction(p + "_up", reversible=False)
    Reac.add_substrates(Subs, stoichiometry=1)
    Reac.add_products(Prod, stoichiometry=1)
    model.add_reactions([Reac])
    model.reactions_info.pop("glk")
    # replacement of the remove_reactions function provided by the package,
    # as it does not work
    return model


# Function to enable CO2 exchange in a single metabolite
def carb_exch(model, s: str, carbs: str, c_unl: str):

    prod_carbs = carbs.replace(c_unl, "")
    Subs = Metabolite(s, atoms=carbs)
    Prod = Metabolite("exchange_partner", atoms=prod_carbs)
    CO2 = Metabolite("CO2", atoms=c_unl)
    Reac = Reaction(s + "_exch", reversible=True)
    Reac.add_substrates(Subs, stoichiometry=1)
    Reac.add_products([Prod, CO2], stoichiometry=[1, 1])
    model.add_reactions([Reac])
    return model


# Real measurements are only given in mdvs (as long as no fragments are measured).
# Hence, we need to consider all possible orientations of labelign patterns within
# the metabolite. We do so by iterating through all possible combinations of
# labeling pattern
def shuffle_labeling_pattern(labeling_strategy: dict):
    permut = [
        set(["".join(j) for j in permutations(l, len(l))])
        for l in labeling_strategy["labeling_pattern"]
    ]
    a = [list(p) for p in product(*permut)]
    labeling_strategies = [
        dict(
            labeled_substrate=labeling_strategy["labeled_substrate"],
            labeling_pattern=x,
            percentage=labeling_strategy["percentage"],
            purity=labeling_strategy["purity"],
        )
        for x in a
    ]
    return labeling_strategies


def get_mdv_results_from_simulations(
    model, fluxes_df, target_emus, labeling_strategy, parallel=None, as_df=True
):
    if parallel is not None:
        vals = [
            (model, fluxes, target_emus, labeling_strategy, sim_name)
            for sim_name, fluxes in fluxes_df.to_dict().items()
        ]
        pool = mp.Pool(parallel)
        all_mdv_res = pool.starmap(get_single_sim_result_from_simulation, vals)
    else:
        all_mdv_res = [
            get_single_sim_result_from_simulation(
                model, fluxes, target_emus, labeling_strategy, sim_name
            )
            for sim_name, fluxes in tqdm(fluxes_df.to_dict().items())
        ]

    if as_df:
        df_freeflux_mdvs = pd.DataFrame(
            [get_mdv_dict(x) for x in all_mdv_res],
            index=[f"s{i:03}" for i in range(1, len(all_mdv_res) + 1)],
        ).transpose()
        return df_freeflux_mdvs
    else:
        return all_mdv_res


def get_single_sim_result_from_simulation(
    model, fluxes, target_emus, labeling_strategy, sim_name="demo"
):
    sim = model.simulator("ss")

    sim.set_target_EMUs(target_emus)
    sim.set_labeling_strategy(**labeling_strategy)
    for r_name, flux in fluxes.items():
        sim.set_flux(r_name, flux)

    sim.prepare()
    res = sim.simulate()
    return res


# function to calculate net_fluxes from fwd and rev
# our current workflow simulates forward and reverse reaction flux as separate values
# for the flux ratio calculation, we have to use the net flux of the reaction
# input: fluxes simulated from sample_fluxes
# returns: dataframe of net fluxes


def calc_net_fluxes(fluxes):
    reactions = fluxes.index
    net_reactions = [re.sub(r"(_f|_b)$", "", r) for r in reactions]
    net_reactions = list(dict.fromkeys(net_reactions))
    net_fluxes = pd.DataFrame(index=net_reactions, columns=fluxes.columns)
    for r in reactions:
        if "_f" in r:
            fwd = fluxes.loc[r]
            rev = fluxes.loc[r.replace("_f", "_b")]
            net = fwd - rev
            net_fluxes.loc[r.replace("_f", "")] = net
        elif "_bm" in r:
            net_fluxes.loc[r] = fluxes.loc[r]
        elif "_b" in r:
            pass
        else:
            net_fluxes.loc[r] = fluxes.loc[r]
    return net_fluxes


# function to build metabolic nodes flux ratios so we know which fluxes we need to use to calculate flux ratios


def build_flux_ratios(net_fluxes, stoich_mat, ratios_of_interest, reactions_to_exclude):
    metabs = net_fluxes.index

    result = {}
    for idx, row in stoich_mat.iterrows():
        pos = row[row == 1].index.tolist()
        pos_renamed = [
            x.replace("_f", "")
            for x in pos
            if "_b" not in x
            if x not in reactions_to_exclude
        ]
        neg = row[row == -1].index.tolist()
        neg_renamed = [
            x.replace("_f", "")
            for x in neg
            if "_b" not in x
            if x not in reactions_to_exclude
        ]
        result[idx] = {"producing": pos_renamed, "consuming": neg_renamed}

    result_filt = {k: result[k] for k in ratios_of_interest}

    return result_filt


# function to calculate flux ratio values in a dataframe

# formulate flux ratios and corresponding formulas


def calc_ratios_from_sims(net_fluxes, flux_ratios):

    ratios = []

    for metab in flux_ratios.keys():
        ratio_output = {}
        n_input = len(flux_ratios[metab]["producing"])
        n_output = len(flux_ratios[metab]["consuming"])
        input_ratios = [
            r + "/" + "+".join(flux_ratios[metab]["producing"])
            for r in flux_ratios[metab]["producing"]
            if n_input > 1
        ]
        output_ratios = [
            r + "/" + "+".join(flux_ratios[metab]["consuming"])
            for r in flux_ratios[metab]["consuming"]
            if n_output > 1
        ]
        ratios += input_ratios
        ratios += output_ratios

    # create an empty dataframe of dimensions flux_ratios x simulations

    flux_ratios_df = pd.DataFrame(index=ratios, columns=net_fluxes.columns)

    for ratio in flux_ratios_df.index:
        num = ratio.split("/")[0]
        denom = ratio.split("/")[1]
        denoms = denom.split("+")
        denom_vals = net_fluxes.loc[denoms]
        num_vals = net_fluxes.loc[num]
        num_vals[num_vals < 0.0001] = 0.0001
        denom_vals[denom_vals < 0.0001] = 0.0001
        ratios = num_vals / denom_vals.sum()
        flux_ratios_df.loc[ratio] = pd.Series(ratios)

    return flux_ratios_df


# Visual inspection of flux ratio values based on flux simulations


def plot_flux_ratios(flux_ratios_df):

    flux_ratios_df["ratio"] = flux_ratios_df.index

    # Melt to long format (so each iteration becomes a row)
    flux_ratios_long = flux_ratios_df.melt(
        id_vars="ratio", var_name="iteration", value_name="value"
    )
    bxplt = sns.catplot(
        data=flux_ratios_long,
        x="ratio",
        y="value",
        kind="box",
        col="ratio",
        sharey=False,
        sharex=False,
        col_wrap=5,
        showfliers=False,
    )

    # Overlay stripplots on each facet
    for ax, ratio_val in zip(bxplt.axes.flatten(), flux_ratios_long["ratio"].unique()):
        ax.set_ylim(0, 1)
        ax.margins(y=0)
        sns.stripplot(
            data=flux_ratios_long[flux_ratios_long["ratio"] == ratio_val],
            x="ratio",
            y="value",
            ax=ax,
            color="black",
            dodge=True,
            jitter=True,
        )
    bxplt.set_titles("")
    bxplt.set_xlabels("")


# given a dictionary of flux ratios with formulas and a dataframe with fluxes, calculate flux ratios and return the calculations in a wide format
def get_flux_ratio_wide(data, nominator, denominator, name):
    nom_fluxes = []
    denom_fluxes = []

    for r_id in nominator:
        if r_id + "_f" in data.index and r_id + "_b" in data.index:
            nom_fluxes.append(data.loc[r_id + "_f"] - data.loc[r_id + "_b"])
        else:
            nom_fluxes.append(data.loc[r_id])

    for r_id in denominator:
        if r_id + "_f" in data.index and r_id + "_b" in data.index:
            denom_fluxes.append(data.loc[r_id + "_f"] - data.loc[r_id + "_b"])
        else:
            denom_fluxes.append(data.loc[r_id])

    ratio = np.sum(nom_fluxes, axis=0) / np.sum(denom_fluxes, axis=0)
    print(ratio)
    return pd.DataFrame([ratio], index=[name], columns=data.columns)


# function to train and validate a flux ratio predictor rf and plot its performance
def eval_RF_sim_pred(
    path_to_fluxes,
    path_to_mdvs,
    path_to_output_dir,
    flux_ratios_dict,
    labeling_strategy: str,
    ratio: str,
):

    fluxes = pd.read_csv(path_to_fluxes, index_col="Unnamed: 0")
    mdvs = pd.read_csv(path_to_mdvs, index_col="Unnamed: 0")

    if list(fluxes.columns) != list(mdvs.columns):
        raise ValueError("Dataframes do not have the same column names!")

    # reimplement so it immediately comes as wide frame

    ratios = pd.concat(
        [
            fluxes.pipe(get_flux_ratio_wide, *reac_ratio, name=name)
            for name, reac_ratio in flux_ratios_dict.items()
        ],
        axis=0,
    )

    # put dfs in correct format for rf
    ratios = ratios.transpose()
    mdvs = mdvs.transpose()

    # select single flux ratio of interest
    ratios = ratios[[ratio]]

    # filter away samples where flux ratio is not between 0 and 1

    # ratios_filt = ratios[ratios[ratio].between(0,1)]
    # mdvs_filt = mdvs.loc[ratios_filt.index]

    print("The number of retained samples is: ", mdvs.shape[0])

    # split data into training and test set

    mdv_train, mdv_test, ratios_train, ratios_test = train_test_split(
        mdvs, ratios, test_size=1 / 3, random_state=78
    )

    # evaluate performance by five-fold cross-validation

    cval = KFold(n_splits=5, shuffle=True, random_state=78)

    rf = RandomForestRegressor(random_state=78)

    scores = cross_val_score(
        rf,
        mdv_train,
        ratios_train.values.ravel(),
        cv=cval,
        scoring="neg_mean_absolute_error",
    )

    score_mean = np.mean(scores)
    score_sd = np.std(scores)

    rf.fit(mdv_train, ratios_train.values.ravel())

    ratios_pred = rf.predict(mdv_test)
    mae = mean_absolute_error(ratios_pred, ratios_test)
    rf.mae_ = mae  # attach mean absolute error to the object to inspect it after export
    pred_sim = pd.DataFrame(
        {"prediction": ratios_pred, "simulation": ratios_test[ratio]}
    )
    fig, ax = plt.subplots(figsize=(6,6))
    ax = sns.scatterplot(data=pred_sim, x="simulation", y="prediction")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(
        0.05, 0.93, "MAE: " + str(round(mae, ndigits=3)), fontsize=22, color="black"
    )
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_aspect('equal', adjustable='box')  
    ax.tick_params(labelsize = 22)
    ax.set_xlabel("Simulation", fontsize = 25)
    ax.set_ylabel("Prediction", fontsize = 25)
    ax.axline((0, 0), slope=1, color='black', linestyle = "--")


    ratio_str = ratio.replace("/", "_")
    plt.savefig(path_to_output_dir + ratio_str + "_" + labeling_strategy + ".pdf")

    plt.show()

    return rf


# subsample simulation iterations per flux ratio to make the training data as even as possible


def sample_ratio_specific_fluxsets(
    path_to_fluxes, path_to_output_dir, flux_ratios_dict, samples_per_bin
):

    # read in fluxes
    fluxes = pd.read_csv(path_to_fluxes, index_col="Unnamed: 0")

    bins = np.linspace(0, 1, 11)

    # calculate flux ratios from fluxes
    flux_ratios = pd.concat(
        [
            pd.concat(
                [
                    fluxes.pipe(get_flux_ratio, *reac_ratio)
                    .assign(flux_ratio_name=name)
                    .reset_index()
                    .rename(columns={"index": "iteration"})
                    for name, reac_ratio in flux_ratios_dict.items()
                ],
                ignore_index=True,
            )
        ]
    )

    # bin flux ratios
    flux_ratios = flux_ratios.assign(
        bin=pd.cut(
            flux_ratios["flux_ratio"], bins=bins, include_lowest=True, right=False
        )
    )
    # subsample a new dataframe which stratifies the training data as good as possible given the input flux distribution
    flux_ratios_balanced = flux_ratios.groupby(
        ["bin", "flux_ratio_name"], group_keys=False
    ).apply(lambda group: sample_if_enough(group, samples_per_bin))

    fluxes.columns = list(range(0, fluxes.shape[1]))

    # for each flux ratio provided to the framework, subset the original flux dataset and save as separate flux distribution file
    for r in flux_ratios_dict.keys():

        flux_ratio_balanced = flux_ratios_balanced.query("flux_ratio_name == @r")

        # get the index from the flux_ratios balanced
        iterations_selected = flux_ratio_balanced["iteration"].values

        # subset the flux dataset for the specific flux ratio
        fluxes_selected = fluxes[iterations_selected]

        fluxes_selected.to_csv(path_to_output_dir + r + "_fluxes.csv")


def get_flux_ratio_onlypos(data, nominator, denominator, name="flux_ratio"):
    nom_fluxes = []
    denom_fluxes = []
    for r_id in nominator:
        if r_id + "_f" in data.index and r_id + "_b" in data.index:
            nom_net = data.loc[r_id + "_f"] - data.loc[r_id + "_b"]
            nom_fluxes.append(nom_net.where(nom_net > 0, np.nan))
        else:
            nom_fluxes.append(data.loc[r_id])
    for r_id in denominator:
        if r_id + "_f" in data.index and r_id + "_b" in data.index:
            denom_net = data.loc[r_id + "_f"] - data.loc[r_id + "_b"]
            denom_fluxes.append(denom_net.where(denom_net > 0, np.nan))
        else:
            denom_fluxes.append(data.loc[r_id])

    return pd.DataFrame(
        np.sum(nom_fluxes, axis=0) / np.sum(denom_fluxes, axis=0), columns=[name]
    )

def PCA_ratiomapping(
    path_to_input_allfluxes: str, 
    path_to_ratio_specific_simulations: str, 
    subsample: int, 
    flux_ratio_dict: dict
):
    
    # calculate principal components
    fluxes = pd.read_csv(path_to_input_allfluxes, index_col = "Unnamed: 0").sample(subsample, axis = 1)
    net_fluxes = calc_net_fluxes(fluxes)
    net_fluxes_scaled = StandardScaler().fit_transform(net_fluxes.transpose())
    pca = PCA(n_components = 2)

    components = pca.fit_transform(net_fluxes_scaled)
    components = pd.DataFrame(components, columns = ["PC1", "PC2"]).reset_index(names = "index")    

    # extract label information for 
    
    all_colnames = []
    for r in flux_ratio_dict.keys():
        fluxes_subs = pd.read_csv(path_to_ratio_specific_simulations + r + "_fluxes.csv", index_col = "Unnamed: 0")
        # get the simulation indices 
        its = pd.DataFrame({
            "ratio": r, 
            "index": fluxes_subs.columns
        })

        all_colnames.append(its)

    all_colnames = pd.concat(all_colnames)
    mask = all_colnames['index'].duplicated(keep=False)  # True for ALL occurrences
    all_colnames.loc[mask, 'ratio'] = 'multiple'
    all_colnames = all_colnames.drop_duplicates()
    all_colnames['ratio'] = all_colnames['ratio'].replace("NaN", "none")


    # bind together 
    all_colnames['index'] = all_colnames['index'].astype(int)
    components_merged = components.merge(all_colnames, on = "index", how = "left")
    components_merged['ratio'] = components_merged['ratio'].fillna("none").astype("category")

    # plot PCA
    
    colors = sns.color_palette('tab10', 14) + [mcolors.to_rgb('gray')]
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data = components_merged, x = "PC1", y = "PC2", hue = "ratio", palette = colors, alpha = 0.5)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    plt.tight_layout()
    plt.savefig(path_to_ratio_specific_simulations + "PCA.pdf", bbox_inches = "tight")
    plt.show()

