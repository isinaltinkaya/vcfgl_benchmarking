###############################################################################
# Step 5: Collect and plot genotype concordance metrics
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

simulation_id = "sim_vcfgl_2310"


configfile: "config/" + simulation_id + ".yaml"


VCFGL = config["tools"]["vcfgl"]
MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]


DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]

BETA_VARS = config["beta_variance_values_neg_e"]

QSERR=['1_betavar'+str(x) for x in BETA_VARS]

BCFTOOLS=config['tools']['bcftools']
PICARDJAR=config['tools']['picardjar']


GTCHECK_ERR=config['bcftools_gtcheck_error']

ploidy=2
DEF_POPS=config['def_pops']

haplo_list=[]
indv_names=[]

for i, (key, value) in enumerate(DEF_POPS.items()):
	haplo_list.append(value*ploidy)
	for ind in range(value):
		# Using PLINK-like format: <Family-ID>_<Individual-ID>
		# to store <Population-ID>_<Individual-ID>
		indv_names.append(f"pop{key}_ind{str(ind+1)}")


###############################################################################
# BEGIN RULES



rule all:
	input:
		# expand(
		# 	"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}_tidy.tsv",
		# 	simid=simulation_id,
		# 	model_id=MODEL,
		# 	contig=CONTIG,
		# 	rep=REP,
		# 	depth=DEPTH,
		# 	error_rate=ERROR_RATE,
		# 	qserr=QSERR,
		# 	gtcheck_err=GTCHECK_ERR,
		# ),
		expand(
			"sim/{simid}/model_{model_id}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}.tsv",
			simid=simulation_id,
			model_id=MODEL,
			gtcheck_err=GTCHECK_ERR,
		)

###############################################################################

# rule collect_picard_get_genotype_concordance:
# 	input:
# 		summary_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_summary_metrics",
# 		contingency_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_contingency_metrics",
# 		genotype_concordance_detail_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_detail_metrics",
# 	output:
# 	shell:
# 		"""
# 		"""




###############################################################################


def get_qs_error(wildcards):
	l=wildcards.qserr.split("1_betavar")
	if(len(l)==2):
		return l[1]
	else:
		return "0"

rule tidy_bcftools_gtcheck:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.txt"
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}_tidy.tsv",
	params:
		qs_error=get_qs_error
	shell:
		"""
		nCompared=$(grep "^INFO" {input}|grep "sites-compared" |cut -f3)
		nNoMatch=$(grep "^INFO" {input}|grep "sites-skipped-no-match" |cut -f3)

		awk -v CONTIG={wildcards.contig} -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v ERRORRATE={wildcards.error_rate} -v QSERR={params.qs_error} -v NCOMPARED=${{nCompared}} -v NNOMATCH=${{nNoMatch}} 'BEGIN{{FS="\t";OFS="\t"}}{{if($1=="DC"){{print $4,$5,$6,CONTIG,REP,DEPTH,ERRORRATE,QSERR,NCOMPARED,NNOMATCH}}}}' {input} > {output}
		"""

rule collect_bcftools_gtcheck:
	input:
		expand(
			"sim/{{simid}}/model_{{model_id}}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{{gtcheck_err}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}_tidy.tsv",
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			qserr=QSERR,
		),
	output:
		"sim/{simid}/model_{model_id}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}.tsv"
	shell:
		"""
		cat {input} > {output}
		"""
