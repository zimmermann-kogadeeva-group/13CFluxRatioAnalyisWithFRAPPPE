from functools import partial

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def xy_line(axes, ls="--", color="black", **kwargs):
    x0, x1 = axes.get_xlim()
    y0, y1 = axes.get_ylim()
    axes.plot(*[[min(x0, y0), max(x1, y1)]] * 2, ls=ls, color=color, **kwargs)


def xy_line_plotly(figure, x0=0, y0=1, x1=1, y1=1, **kwargs):
    if "line" not in kwargs:
        kwargs["line"] = dict(dash="dash")
    if "row" not in kwargs:
        kwargs["row"] = "all"
    if "col" not in kwargs:
        kwargs["col"] = "all"
    figure.add_shape(type="line", x0=x0, y0=y0, x1=x1, y1=y1, **kwargs)


def melt_flux_df(data, **kwargs):
    return (
        data.reset_index(names="reaction_id")
        .melt(id_vars="reaction_id", var_name="sim_name", value_name="flux")
        .assign(**kwargs)
    )


def melt_mdv_df(data, **kwargs):
    return (
        data.melt(ignore_index=False, var_name="sim_name", value_name="mdv")
        .reset_index(names="emu")
        .pipe(
            lambda x: x.join(x.emu.str.extract(r"(?P<mb>.*)_(?P<atomnos>.*)_(?P<m>.*)"))
        )
        .assign(emu=lambda x: x.mb + "_" + x.atomnos)
    )


def get_sumoflux_fluxes(file_prefix):
    df_sumoflux_f = pd.read_csv(f"{file_prefix}_f.csv")
    df_sumoflux_b = pd.read_csv(f"{file_prefix}_b.csv")
    return pd.concat(
        [
            df_sumoflux_f.query("r_irreversible == 1"),
            df_sumoflux_b.query("r_irreversible == 0").assign(
                Row=lambda x: x.Row + "_b"
            ),
            df_sumoflux_f.query("r_irreversible == 0").assign(
                Row=lambda x: x.Row + "_f"
            ),
        ],
        ignore_index=True,
    )


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


def get_highest_diff(var_name, value_name, hopsy_col, sumoflux_col, n=10):
    return (
        (sumoflux_col - hopsy_col)
        .sort_values(ascending=False, key=np.absolute)
        .head(n)
        .reset_index(name=value_name)
        .rename(columns={"index": var_name})
    )


get_highest_flux_diff = partial(get_highest_diff, "reaction_id", "flux_diff")
get_highest_mdv_diff = partial(get_highest_diff, "emu", "mdv_diff")


def convert_sumoflux_notation(metabolite):
    met, suf = metabolite.split("_", maxsplit=2)
    min_val = 1
    if suf == met:
        suf = 0

    if met[-2].isdigit():
        min_val, max_val = int(met[-2]), int(met[-1])
        new_labelling = "".join([str(x) for x in range(int(min_val), int(max_val) + 1)])
        name = met[:-2]
        if suf != 0:
            suf = int(suf) - min_val + 1
    else:
        new_labelling = met[-1]
        name = met[:-1]
        if suf != 0:
            suf = 1

    return f"{name}_{new_labelling}_{suf}"


def plot_fluxes(func=None, **kwargs):
    figsize = (10, 16) if "figsize" not in kwargs else kwargs.pop("figsize")
    func = func or sns.stripplot

    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(right=0.80)
    func(
        pd.concat(
            [
                ind_flux_df.pipe(melt_flux_df, comparison=name)
                for name, ind_flux_df in kwargs.items()
            ]
        ),
        y="reaction_id",
        x="flux",
        hue="comparison",
        dodge=True,
    )
    sns.move_legend(ax, loc="center left", bbox_to_anchor=(1.0, 0.5))
    return fig


def plot_mdvs_boxplots(data, col_order=None, **kwargs):
    fig = sns.FacetGrid(
        data.pipe(melt_mdv_df), col="emu", col_wrap=4, col_order=col_order
    ).map_dataframe(sns.boxplot, x="m", y="mdv", **kwargs)
    return fig


def plot_rf_results(
    results_dict, col_wrap=2, sharex=False, sharey=False, height=5, **kwargs
):
    df_pred_sim = pd.concat(
        [
            results[2].assign(flux_ratio_name=fr_name)
            for fr_name, results in results_dict.items()
        ],
        ignore_index=True,
    )

    fg = sns.FacetGrid(
        df_pred_sim,
        col="flux_ratio_name",
        col_wrap=col_wrap,
        sharex=sharex,
        sharey=sharey,
        height=height,
        **kwargs,
    ).map_dataframe(
        sns.scatterplot,
        x="simulation",
        y="prediction",
    )
    for fr_name, ax in fg.axes_dict.items():
        scores = results_dict[fr_name][1]
        ax.text(
            0.05,
            0.93,
            f"MAE: {scores['mae']:.3f}",
            fontsize=10,
            color="black",
            transform=ax.transAxes,
        )
        ax.text(
            0.05,
            0.88,
            f"CV MAE mean: {scores['cv_mae_mean']:.3f} \u00b1 {scores['cv_mae_sd']:.3f}",
            fontsize=10,
            color="black",
            transform=ax.transAxes,
        )
        xy_line(ax)

    return fg
