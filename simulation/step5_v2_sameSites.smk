###############################################################################
# Step 5: Collect and plot genotype discordance metrics
#
# isinaltinkaya
###############################################################################
#
# Naming convention:
# - log files:
# output: sim/{simid}/REST_OF_PATH.EXTENSION
# params.prefix: sim/{simid}/REST_OF_PATH
# log: sim/{simid}/logs/REST_OF_PATH.EXTENSION

import stdpopsim, msprime, tskit
import numpy as np

simulation_id = config["simulation_id"]



MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]

DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]

ploidy = 2
DEF_POPS = config["def_pops"]

haplo_list = []
indv_names = []

for i, (key, value) in enumerate(DEF_POPS.items()):
	haplo_list.append(value * ploidy)
	for ind in range(value):
		# Using PLINK-like format: <Family-ID>_<Individual-ID>
		# to store <Population-ID>_<Individual-ID>
		indv_names.append(f"pop{key}_ind{str(ind+1)}")


###############################################################################
# BEGIN RULES


GC_METHODS = ["genotype_calling", "genotype_calling_perpop"]


rule all:
	input:
		expand(
			"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ{dogq}.sameSites.tsv.tidy",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			gc_method=GC_METHODS,
			dogq=[5, 6],
			glSpecs=["precise1_gl2", "gl1", "gl2"],
			qsbeta=["0_0", "2_5", "2_6"],
		),
		expand(
			"sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}-doGQ{dogq}_sameSites.tsv",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			gc_method=GC_METHODS,
			dogq=[5, 6],
			glSpecs=["precise1_gl2", "gl1", "gl2"],
		),


rule tidy_get_genotype_discordance_doGQ:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}-snp_doGQ{dogq}.sameSites.tsv",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}-snp_doGQ{dogq}.sameSites.tsv.tidy",
	shell:
		"""
		awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v BETAVAL={wildcards.betaval} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,BETAVAL}}' {input} > {output}
		"""


rule collect_get_genotype_discordance_qs_doGQ:
	input:
		expand(
			"sim/{{simid}}/model_{{model_id}}/contig_{contig}/gc_evaluation/genotype_discordance/{{gc_method}}_{{glSpecs}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ{{dogq}}.sameSites.tsv.tidy",
			contig=CONTIG,
			depth=DEPTH,
			rep=REP,
			error_rate=ERROR_RATE,
			qsbeta=["0_0", "2_5", "2_6"],
		),
	output:
		"sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}-doGQ{dogq}_sameSites.tsv",
	shell:
		"""
		cat {input} > {output}
		"""

