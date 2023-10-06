###############################################################################
# Step 4: Evaluate genotype calling
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

BCFTOOLS=config['tools']['bcftools']
PICARDJAR=config['tools']['picardjar']


MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]
DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]
BETA_VARS = config["beta_variance_values_neg_e"]
QSERR=['1_betavar'+str(x) for x in BETA_VARS]

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
# # # #1
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),
# # #9
		expand("sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/get_genotype_discordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.tsv",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			qserr=QSERR,
		),
# # # # #2
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),
# # # # #3
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),
# # # # #4
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.txt",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 			gtcheck_err=GTCHECK_ERR,
# 		),




# ###############################################################################

# # # #1
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),

# TODO validate that the true genotypes are all the same across different qserrs
rule get_real_gt:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
	shell:
		"""
		{BCFTOOLS} annotate -x "FORMAT/GL,FORMAT/PL,FORMAT/GP,FORMAT/DP,FORMAT/AD,INFO/QS,INFO/AD,INFO/AC,INFO/AN" {input} |{BCFTOOLS} view --trim-alt-alleles --no-update -Ob -o {output}
		"""


# # ###############################################################################

# # # # # #2
# # 		expand(
# # 			"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
# # 			simid=simulation_id,
# # 			model_id=MODEL,
# # 			contig=CONTIG,
# # 			rep=REP,
# # 			depth=DEPTH,
# # 			error_rate=ERROR_RATE,
# # 			qserr=QSERR,
# # 		),


# rule index_for_bcftools_gtcheck_callgt:
# 	input:
# 		"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 	output:
# 		"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
# 	shell:
# 		"""
# 		{BCFTOOLS} index {input};
# 		"""

# # ###############################################################################

# # # # # #3
# # 		expand(
# # 			"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
# # 			simid=simulation_id,
# # 			model_id=MODEL,
# # 			contig=CONTIG,
# # 			rep=REP,
# # 			depth=DEPTH,
# # 			error_rate=ERROR_RATE,
# # 			qserr=QSERR,
# # 		),

# rule index_for_bcftools_gtcheck_truthgt:
# 	input:
# 		"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 	output:
# 		"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
# 	shell:
# 		"""
# 		{BCFTOOLS} index {input};
# 		"""


# # ###############################################################################

# # # # # #4
# # 		expand(
# # 			"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.txt",
# # 			simid=simulation_id,
# # 			model_id=MODEL,
# # 			contig=CONTIG,
# # 			rep=REP,
# # 			depth=DEPTH,
# # 			error_rate=ERROR_RATE,
# # 			qserr=QSERR,
# # 			gtcheck_err=GTCHECK_ERR,
# # 		),

# rule bcftools_gtcheck_noerror:
# 	input:
# 		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 		callgtcsi="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
# 		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 		truthgtcsi="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
# 	output:
# 		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.txt",
# 	params:
# 		pvar=str(','.join([i+","+i for i in indv_names]))
# 	log:
# 		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.txt",
# 	shell:
# 		"""
# 		({BCFTOOLS} gtcheck -u GT,GT -p {params.pvar} -g {input.truthgt} {input.callgt} -e {wildcards.gtcheck_err} > {output})2>{log}
# 		"""


# ###############################################################################

# # # #5
		# expand(
		# 	"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps-sample{sid}.vcf",
		# 	simid=simulation_id,
		# 	model_id=MODEL,
		# 	contig=CONTIG,
		# 	rep=REP,
		# 	depth=DEPTH,
		# 	error_rate=ERROR_RATE,
		# 	qserr=QSERR,
		# 	sid=indv_names,
		# ),

# rule bcf_to_indvcf_for_picard_callgt:
# 	input:
# 		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf",
# 	output:
# 		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps-sample{sid}.vcf",
# 	shell:
# 		"""
# 		{BCFTOOLS} view -O v -s {wildcards.sid} -o {output.callgt} {input.callgt}
# 		"""


# # # #6
# 		expand(
# 	 		"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.vcf",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 			sid=indv_names,
# 		),

# rule bcf_to_indvcf_for_picard_truthgt:
# 	input:
# 		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 	output:
# 		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.vcf",
# 	shell:
# 		"""
# 		{BCFTOOLS} view -O v -s {wildcards.sid} -o {output.truthgt} {input.truthgt}
# 		"""



# # # #7
# 		expand(
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 			sid=indv_names,
# 		),

# rule picard_get_genotype_concordance:
# 	input:
# 		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.vcf",
# 		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps-sample{sid}.vcf",
# 	output:
# 		summary_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.genotype_concordance_summary_metrics",
# 		contingency_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.genotype_concordance_contingency_metrics",
# 		genotype_concordance_detail_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.genotype_concordance_detail_metrics",
# 	params:
# 		prefix="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}",
# 	log:
# 		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.genotype_concordance_summary_metrics",
# 	shell:
# 		"""
# 		java -jar {PICARDJAR} GenotypeConcordance \
# 				CALL_VCF={input.callgt} \
# 				CALL_SAMPLE={wildcards.sid} \
# 				O={params.prefix} \
# 				TRUTH_VCF={input.truthgt} \
# 				TRUTH_SAMPLE={wildcards.sid} \
# 				VERBOSITY=WARNING
# 		done
# 		"""


# # # #8
# 		expand(
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 			sid=indv_names,
# 		),

# rule merge_persample_picard_genotype_concordance:
# 	input:
# 		summary_metrics=expand("sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.genotype_concordance_summary_metrics", sid=indv_names),
# 		contingency_metrics=expand("sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.genotype_concordance_contingency_metrics", sid=indv_names),
# 		genotype_concordance_detail_metrics=expand("sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-sample{sid}.genotype_concordance_detail_metrics", sid=indv_names),
# 	output:
# 		summary_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_summary_metrics",
# 		contingency_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_contingency_metrics",
# 		genotype_concordance_detail_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_detail_metrics",
# 	shell:
# 		"""
# 		cat {input.summary_metrics} > {output.summary_metrics}
# 		cat {input.contingency_metrics} > {output.contingency_metrics}
# 		cat {input.genotype_concordance_detail_metrics} > {output.genotype_concordance_detail_metrics}
# 		"""

###############################################################################

# # # #9
# 		expand(
# 		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/get_genotype_discordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.tsv",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),

GET_GT_DISCORDANCE="/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/simulation/misc/get_genotype_discordance"

rule get_genotype_discordance:
	input:
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/get_genotype_discordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.tsv",
	params:
		pvar=str(','.join([i+","+i for i in indv_names])),
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/get_genotype_discordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.tsv",
	shell:
		"""
		({GET_GT_DISCORDANCE} {input.truthgt} {input.callgt} > {output})2> {log}
		"""
