###############################################################################
# Step 3: Run vcfgl simulation with qs error
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

simulation_id= "sim_vcfgl_2309"
configfile: "config/"+simulation_id+".yaml"

MODEL=config['model']
ERROR_RATE=config['vcfgl_error_rate']
DEPTH=config['depth']
CONTIG=config['contig']
REP = [*range(config['n_reps'])]

BETA_VARS=config['beta_variance_values_neg_e']
QSERR=['1_betavar'+str(x) for x in BETA_VARS]

VCFGL=config['tools']['vcfgl']

###############################################################################
# BEGIN RULES

rule all:
	input:
		expand("sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf",
				simid=simulation_id,
				model_id=MODEL,
				contig=CONTIG,
				rep=REP,
				depth=DEPTH,
				error_rate=ERROR_RATE,
				qserr=QSERR,
				),



rule simulate_vcfgl_qs_error:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
	params:
		prefix="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}",
		random_seed=lambda wildcards: str(int(wildcards.rep)+1),
		error_qs=lambda wildcards: str(wildcards.qserr).split('_')[0],
		beta_var=lambda wildcards: str(wildcards.qserr).split('_')[1][7],
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}.bcf"
	shell:
		"""
		({VCFGL} -i {input} \
				-o {params.prefix} \
				-O b \
				--depth {wildcards.depth} \
				--error-rate {wildcards.error_rate}  \
				--error-qs {params.error_qs} \
				--beta-variance 1e-{params.beta_var} \
				-addQS 1 \
				-explode 1 \
				-addFormatAD 1 \
				-addInfoAD 1 \
				-addPL 1 \
				-addGP 1 \
				--seed {params.random_seed} ) 2> {log}
		"""

