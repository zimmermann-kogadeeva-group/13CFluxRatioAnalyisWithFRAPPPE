
import pandas as pd
from pathlib import Path

df_runs = pd.read_csv("Config/runsheet.tsv", sep="\t", index_col=0)


rule all:
	input:
		expand("{output_dir}/fluxes.csv", output_dir=df_runs.index)


def get_sim_inputs(wildcards):
	outdir = f"Output/{wildcards.path_to_run}"
	return df_runs.loc[outdir, :].to_dict()
	

rule runs:
	input:
		unpack(get_sim_inputs)
	output:
		"Output/{path_to_run}/fluxes.csv",
		"Output/{path_to_run}/mdvs.csv"
	conda:
		"Envs/runs.yaml"
	params:
		output_dir = lambda wc, output: str(Path(output[0]).parent),
		script = lambda wc, input: input.script,
		parallel = 25
	resources:
		cpus = 25,
		mem_mb = lambda wc, attempt: 2400 * attempt,
		runtime = lambda wc, attempt: 240 * attempt
	script:
		"{params.script}"

