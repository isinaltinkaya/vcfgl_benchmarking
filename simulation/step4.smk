###############################################################################
# Step 3: SNP calling and genotype calling
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

MODEL=config['model']
ERROR_RATE=config['vcfgl_error_rate']
DEPTH=config['depth']
#TODO 
# DEPTH=[0.1]
CONTIG=config['contig']
REP = [*range(config['n_reps'])]
#TODO
REP=0

BETA_VARS=config['beta_variance_values_neg_e']
QSERR=['1_betavar'+str(x) for x in BETA_VARS]

ANGSD=config['tools']['angsd']
BCFTOOLS=config['tools']['bcftools']



###############################################################################
# BEGIN RULES




# TODO DELME

print(config['repset'])
repset=config['repset']
if(repset=="first"):
	REP = [0]
elif(repset=="rest"):
	REP = [*range(1,config['n_reps'])]
elif(repset=="singletest"):
	# sim/sim_vcfgl_2309/model_OutOfAfrica_3G09/contig_chr22/vcfgl/sim_vcfgl_2309-OutOfAfrica_3G09-chr22-rep12-d10-e0.002-qs1_betavar7.bcf
	REP=0
	DEPTH=10
	QSERR='1_betavar7'
else:
	exit("Please specify repset in config file")

print(REP)

rule all:
	input:
# #6
		expand(
			"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			qserr=QSERR,
		),
# # #5
		# expand(
			# "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
			# simid=simulation_id,
			# model_id=MODEL,
			# contig=CONTIG,
			# rep=REP,
			# depth=DEPTH,
			# error_rate=ERROR_RATE,
			# qserr=QSERR,
		# ),
# # #4
		# expand(
			# "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
			# simid=simulation_id,
			# model_id=MODEL,
			# contig=CONTIG,
			# rep=REP,
			# depth=DEPTH,
			# error_rate=ERROR_RATE,
			# qserr=QSERR,
		# ),
		#
# # # #3
# # 		expand(
# # 			"sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.snp_sites.tsv",
# # 			simid=simulation_id,
# # 			model_id=MODEL,
# # 			contig=CONTIG,
# # 			rep=REP,
# # 			depth=DEPTH,
# # 			error_rate=ERROR_RATE,
# # 			qserr=QSERR,
# # 		),
# # # #2
# # 		expand(
# # 			"sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.mafs.gz",
# # 			simid=simulation_id,
# # 			model_id=MODEL,
# # 			contig=CONTIG,
# # 			rep=REP,
# # 			depth=DEPTH,
# # 			error_rate=ERROR_RATE,
# # 			qserr=QSERR,
# # 		),
# # # # #1
# # # 		expand(
# # # 			"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# # # 			simid=simulation_id,
# # # 			model_id=MODEL,
# # # 			contig=CONTIG,
# # # 			rep=REP,
# # # 			depth=DEPTH,
# # # 			error_rate=ERROR_RATE,
# # # 			qserr=QSERR,
# # # 		),


###############################################################################

# # # #1
# # 		expand(
# # 			"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# # 			simid=simulation_id,
# # 			model_id=MODEL,
# # 			contig=CONTIG,
# # 			rep=REP,
# # 			depth=DEPTH,
# # 			error_rate=ERROR_RATE,
# # 			qserr=QSERR,
# # 		),

# rule remove_gts_from_vcfgl_vcfs:
# 	input:
# 		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
# 	output:
# 		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
# 	log:
# 		"sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
# 	shell:
# 		"""
# 		({BCFTOOLS} annotate -x FORMAT/GT {input} -O b -o {output} )2> {log}
# 		"""

###############################################################################

# # #2
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.mafs.gz",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),

rule call_snps_angsd_vcfgl:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.mafs.gz"
	params:
		snppval=config["SNP_calling_SNP_pval"],
		prefix="sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}"
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
	shell:
		"""
		( {ANGSD} \
				-vcf-gl {input} \
				-out {params.prefix} \
				-doMajorMinor 1 \
				-doMaf 1 \
				-skipTriallelic 1 \
				-SNP_pval {params.snppval} \
				-doBcf 1 \
				-doGlf 2 ) 2> {log}
		"""


###############################################################################

# # #3
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.snp_sites.tsv",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),

rule get_snp_calling_snps_bed:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.mafs.gz",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.snp_sites.tsv",
	shell:
		"""
		zcat {input} | sed 1d|cut -f-2 > {output}
		"""


###############################################################################

# # #4
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),

rule bcftools_call_genotypes:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	shell:
		"""
		({BCFTOOLS} call -m -O b -o {output} {input} )2> {log}
		"""


###############################################################################

# # #5
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

rule index_calledgt_bcf:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
	shell:
		"""
		{BCFTOOLS} index {input}
		"""


###############################################################################

# # #6
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),

# angsd removes QS tag, so use angsd sites output to manually filter for snps
rule get_called_genotypes_called_snp_regions:
	input:
		snp_sites="sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.snp_sites.tsv",
		calledbcf="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
		calledbcfcsi="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf.csi",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snps.bcf",
	shell:
		"""
		({BCFTOOLS} view {input.calledbcf} -R {input.snp_sites} -O b -o {output} ) 2> {log}
		"""


