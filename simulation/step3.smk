###############################################################################
# Step 3: Genotype calling
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

MODEL=config['model']
ERROR_RATE=config['vcfgl_error_rate']
DEPTH=config['depth']
CONTIG=config['contig']
REP = [*range(config['n_reps'])]

BETA_VARS=config['beta_variance_values_neg_e']
QSERR=['1_betavar'+str(x) for x in BETA_VARS]

ANGSD=config['tools']['angsd']
BCFTOOLS=config['tools']['bcftools']



###############################################################################
# BEGIN RULES

# print(config['repset'])
# repset=config['repset']
# if(repset=="first"):
# 	REP = [0]
# elif(repset=="rest"):
# 	REP.remove(0)


rule all:
	input:
# # #1
# 		expand(
# 			"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
# 			simid=simulation_id,
# 			model_id=MODEL,
# 			contig=CONTIG,
# 			rep=REP,
# 			depth=DEPTH,
# 			error_rate=ERROR_RATE,
# 			qserr=QSERR,
# 		),
# #2
		expand(
			"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			qserr=QSERR,
		),

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

rule remove_gts_from_vcfgl_vcfs:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
	shell:
		"""
		({BCFTOOLS} annotate -x FORMAT/GT {input} -O b -o {output} )2> {log}
		"""

###############################################################################
# # #2
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
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
	shell:
		"""
		({BCFTOOLS} call -m -O b -o {output} {input} )2> {log}
		"""
