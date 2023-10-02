###############################################################################
# Step 5: Evaluate genotype calling
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

simulation_id = "sim_vcfgl_2309"


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

# TODO DELME

print(config['repset'])
repset=config['repset']
if(repset=="first"):
	REP = [0]
	DEPTH.remove(10)
	#TODO!!! rerun for 10
elif(repset=="rest"):
	REP = [*range(1,config['n_reps'])]
elif(repset=="singletest"):
	REP=0
	DEPTH=10
	# QSERR='1_betavar7'
	QSERR=['1_betavar7','1_betavar6']

else:
	exit("Please specify repset in config file")

print(REP)


rule all:
	input:
		# # #3
		expand(
			"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.txt",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			qserr=QSERR,
		),
# #2
	# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_summary_metrics",
	# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_contingency_metrics",
		expand(
			"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_detail_metrics",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			qserr=QSERR,
		),
# #1
		expand(
			"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			qserr=QSERR,
			),
#
#
###############################################################################

# # # #1

# 		expand(
#       "sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),

rule get_real_gt:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
	shell:
		"""
		bcftools annotate -x "FORMAT/GL,FORMAT/PL" {input} |bcftools view --trim-alt-alleles -Ob -o {output}
		"""



###############################################################################

# # # #2

# 		expand(
#     # "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_summary_metrics",
#     # "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_contingency_metrics",
#      "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_detail_metrics",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),

rule bcf_to_vcf_for_picard:
	input:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf",
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	output:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.vcf",
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.vcf",
	shell:
		"""
		{BCFTOOLS} view -O v -o {output.callgt} {input.callgt}
		{BCFTOOLS} view -O v -o {output.truthgt} {input.truthgt}
		"""

rule picard_get_genotype_concordance:
	input:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.vcf",
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.vcf",
	output:
		summary_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_summary_metrics",
		contingency_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_contingency_metrics",
		genotype_concordance_detail_metrics="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_detail_metrics",
	params:
		prefix="sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}",
		samples=indv_names
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/picard_genotype_concordance/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.genotype_concordance_summary_metrics",
	shell:
		"""
		for sid in {params.samples};do
			java -jar {PICARDJAR} GenotypeConcordance \
					CALL_VCF={input.callgt} \
					CALL_SAMPLE=${{sid}} \
					O={output.summary_metrics}_${{sid}}_concordance \
					TRUTH_VCF={input.truthgt} \
					TRUTH_SAMPLE=${{sid}} \
					VERBOSITY=WARNING
			cat {output.summary_metrics}_${{sid}}_concordance.genotype_concordance_summary_metrics >> {output.summary_metrics}
			cat {output.summary_metrics}_${{sid}}_concordance.genotype_concordance_contingency_metrics >> {output.contingency_metrics}
			cat {output.summary_metrics}_${{sid}}_concordance.genotype_concordance_detail_metrics >> {output.genotype_concordance_detail_metrics}
		done
		"""




###############################################################################

# # # #3

# 		expand(
#       "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.txt",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),


rule index_for_bcftools_gtcheck:
	input:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf",
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	output:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf.csi",
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
	shell:
		"""
		{BCFTOOLS} index {input.callgt};
		{BCFTOOLS} index {input.truthgt};
		"""

rule bcftools_gtcheck:
	input:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf",
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
		callgtcsi="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf.csi",
		truthgtcsi="sim/{simid}/model_{model_id}/contig_{contig}/true_genotypes/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.txt"
	params:
		pvar=str(','.join([i+","+i for i in indv_names]))
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/bcftools_gtcheck/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.txt"
	shell:
		"""
		({BCFTOOLS} gtcheck -u GT -p {params.pvar} -g {input.truthgt} {input.callgt} > {output})2>{log}
		"""

###############################################################################
