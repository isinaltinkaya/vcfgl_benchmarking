###############################################################################
# Step 2: Run vcfgl simulation
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

VCFGL=config['tools']['vcfgl']

MODEL=config['model']

ERROR_RATE=config['vcfgl_error_rate']

VCFGL=config['tools']['vcfgl']

DEPTH=config['depth']
CONTIG=config['contig']
REP = [*range(config['n_reps'])]


BETA_VARS=config['beta_variance_values_neg_e']

###############################################################################
# BEGIN RULES

rule all:
	input:
		expand("sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}.bcf",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE
			),
		expand("sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}.bcf",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			beta_variance=BETA_VARS
			)



# SIMULATE GENOTYPE LIKELIHOODS
# no qs error
rule simulate_vcfgl_no_qs_error:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}.bcf"
	params:
		random_seed=lambda wildcards: int(int(wildcards.rep)+1),
		prefix="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}"
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}.bcf"
	shell:
		"""
		({VCFGL} -i {input} \
				-o {params.prefix} \
				-O b \
				--depth {wildcards.depth} \
				--error-rate {wildcards.error_rate}  \
				--error-qs 0 \
				-addQS 1 \
				-explode 1 \
				-addFormatAD 1 \
				-addInfoAD 1 \
				-addPL 1 \
				-addGP 1 \
				--seed {params.random_seed} ) 2> {log}
		"""

rule simulate_vcfgl_qs_error:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}.bcf"
	params:
		random_seed=lambda wildcards: int(int(wildcards.rep)+1),
		prefix="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}"
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}.bcf"
	shell:
		"""
		({VCFGL} -i {input} \
				-o {params.prefix} \
				-O b \
				--depth {wildcards.depth} \
				--error-rate {wildcards.error_rate}  \
				--error-qs 1 \
				--beta-variance 1e-{wildcards.beta_variance} \
				-addQS 1 \
				-explode 1 \
				-addFormatAD 1 \
				-addInfoAD 1 \
				-addPL 1 \
				-addGP 1 \
				--seed {params.random_seed} ) 2> {log}
		"""
