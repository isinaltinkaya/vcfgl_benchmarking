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

# import stdpopsim, msprime, tskit
# import numpy as np

simulation_id = "sim_vcfgl_2312"


configfile: "config/" + simulation_id + ".yaml"


BCFTOOLS = config["tools"]["bcftools"]
PICARDJAR = config["tools"]["picardjar"]


MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]
DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]
BETA_VARS = config["beta_variance_values_neg_e"]

GTCHECK_ERR = config["bcftools_gtcheck_error"]

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


GC_METHODS = ["genotype_calling", "genotype_calling_perpop"]

GQ_FILTERS = ["-minGQ20", ""]
MINGQ = 20

GET_GT_DISCORDANCE = config["tools"]["gtDiscordance"]

GLMODELS = [1, 2]

###############################################################################
# BEGIN RULES

###############################################################################
###############################################################################
rule get_vcf_stats:
	 input:
		 "sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs1_betavar{betavar}-snp.bcf",
	 output:
		 "sim/{simid}/model_{model_id}/contig_{contig}/vcf_stats/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs1_betavar{betavar}-snp_vcf_stats.txt",
	 shell:
		 """
		  {BCFTOOLS} stats {input} > {output}
		  """

rule filter_GQ_for_get_genotype_discordance:
	 input:
		 "sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs1_betavar{betavar}-snp.bcf",
	 output:
		 "sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs1_betavar{betavar}-snp-minGQ{mingq}.bcf",
	 shell:
		 """
		  {BCFTOOLS} +setGT {input} -- -t q -n . -i 'GQ<20' | {BCFTOOLS} annotate -x "INFO/AD,INFO/AC,INFO/AN"| {BCFTOOLS} view --exclude-uncalled -v snps -O b -o {output}
		  """

rule get_genotype_discordance_mingq:
	 input:
		 truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs1_betavar{betavar}.bcf",
		 callgt="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs1_betavar{betavar}-snp-minGQ{mingq}.bcf",
	 output:
		 stats="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs1_betavar{betavar}-snp-minGQ{mingq}_stats.tsv",
		 qs="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs1_betavar{betavar}-snp-minGQ{mingq}_qs.tsv",
	 log:
		 "sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs1_betavar{betavar}-snp-minGQ{mingq}_stats.tsv",
	 shell:
		 """
		 ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output.stats} -doGQ 1 > {output.qs} )2> {log}
		 """

###############################################################################
##############################################################################
## OLD
# ###############################################################################

rule index_for_bcftools_gtcheck_callgt:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf.csi",
	shell:
		"""
		{BCFTOOLS} index {input};
		"""


# ###############################################################################
# # # # #3
rule index_for_bcftools_gtcheck_truthgt:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf.csi",
	shell:
		"""
		{BCFTOOLS} index {input};
		"""

# ###############################################################################
# # # # #4
rule bcftools_gtcheck_noerror:
	input:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
		callgtcsi="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf.csi",
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
		truthgtcsi="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf.csi",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.txt",
	params:
		pvar=str(','.join([i+","+i for i in indv_names]))
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck_e{gtcheck_err}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.txt",
	shell:
		"""
		({BCFTOOLS} gtcheck -u GT,GT -p {params.pvar} -g {input.truthgt} {input.callgt} -e {wildcards.gtcheck_err} > {output})2>{log}
		"""

###############################################################################
# # #5
rule bcf_to_indvcf_for_picard_callgt:
	input:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snps.bcf",
	output:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snps-sample{sid}.vcf",
	shell:
		"""
		{BCFTOOLS} view -O v -s {wildcards.sid} -o {output.callgt} {input.callgt}
		"""

rule bcf_to_indvcf_for_picard_truthgt:
	input:
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
	output:
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}.vcf",
	shell:
		"""
		{BCFTOOLS} view -O v -s {wildcards.sid} -o {output.truthgt} {input.truthgt}
		"""

rule picard_get_genotype_concordance:
	input:
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}.vcf",
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snps-sample{sid}.vcf",
	output:
		summary_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}.genotype_concordance_summary_metrics",
		contingency_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}.genotype_concordance_contingency_metrics",
		genotype_concordance_detail_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}.genotype_concordance_detail_metrics",
	params:
		prefix="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}.genotype_concordance_summary_metrics",
	shell:
		"""
		java -jar {PICARDJAR} GenotypeConcordance \
						CALL_VCF={input.callgt} \
						CALL_SAMPLE={wildcards.sid} \
						O={params.prefix} \
						TRUTH_VCF={input.truthgt} \
						TRUTH_SAMPLE={wildcards.sid} \
						VERBOSITY=WARNING
		done
		"""

rule merge_persample_picard_genotype_concordance:
	input:
		summary_metrics=expand("sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}.genotype_concordance_summary_metrics", sid=indv_names),
		contingency_metrics=expand("sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}.genotype_concordance_contingency_metrics", sid=indv_names),
		genotype_concordance_detail_metrics=expand("sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-sample{sid}.genotype_concordance_detail_metrics", sid=indv_names),
	output:
		summary_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-d{depth}-e{error_rate}-qs{qsbeta}.genotype_concordance_summary_metrics",
		contingency_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-d{depth}-e{error_rate}-qs{qsbeta}.genotype_concordance_contingency_metrics",
		genotype_concordance_detail_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-d{depth}-e{error_rate}-qs{qsbeta}.genotype_concordance_detail_metrics",
	shell:
		"""
		cat {input.summary_metrics} > {output.summary_metrics}
		cat {input.contingency_metrics} > {output.contingency_metrics}
		cat {input.genotype_concordance_detail_metrics} > {output.genotype_concordance_detail_metrics}
		"""
