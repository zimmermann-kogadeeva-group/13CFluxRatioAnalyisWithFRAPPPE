#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

from mfa import (
    create_add_constraints,
    get_mdv_results_from_simulations,
    get_model,
    make_bounds,
    sample_fluxes,
)

try:
    from snakemake.script import snakemake
except ImportError:
    snakemake = None


def read_config(path_to_file):
    with open(path_to_file) as fh:
        return json.load(fh)


def main(
    model_path, bounds, target_emus, labelling, output_dir, flux_ratio_bounds, parallel=None, **kwargs
):
    bounds_dict = read_config(bounds)
    target_emus_dict = read_config(target_emus)
    labelling_dict = read_config(labelling)
    fr_bounds_dict = read_config(flux_ratio_bounds)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model = get_model(model_path)

    stoich_mat = model.get_total_stoichiometric_matrix()

    add_constraints = create_add_constraints(stoich_mat, fr_bounds_dict)

    lb, ub = make_bounds(stoich_mat, **bounds_dict)
    df_fluxes = sample_fluxes(stoich_mat, 50, lb, ub, thinning=10000, add_constraints=add_constraints)
    df_fluxes.to_csv(output_dir / "fluxes.csv")

    df_mdvs = get_mdv_results_from_simulations(
        model,
        df_fluxes,
        target_emus_dict,
        labelling_dict,
        parallel=parallel,
    )
    df_mdvs.to_csv(output_dir / "mdvs.csv")


if __name__ == "__main__":
    if snakemake is not None:
        main(
            **snakemake.input,
            output_dir=snakemake.params["output_dir"],
            parallel=snakemake.threads,
        )
    else:
        parser = argparse.ArgumentParser()
        parser.add_argument("model_path")
        parser.add_argument("bounds")
        parser.add_argument("target_emus")
        parser.add_argument("labelling")
        parser.add_argument("flux_ratio_bounds")
        parser.add_argument("output_dir")
        parser.add_argument("-p", "--parallel", default=1, type=int)
        args = parser.parse_args()
        if args.parallel < 2:
            args.parallel = None

        main(**vars(args))
