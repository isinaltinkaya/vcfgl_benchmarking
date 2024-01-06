###############################################################################
# Step 4: Evaluate genotype calling v2 use same sites
#
# isinaltinkaya
###############################################################################
#
# Naming convention:
# - log files:
# output: sim/{simid}/REST_OF_PATH.EXTENSION
# params.prefix: sim/{simid}/REST_OF_PATH
# log: sim/{simid}/logs/REST_OF_PATH.EXTENSION

simulation_id = config["simulation_id"]



BCFTOOLS = config["tools"]["bcftools"]

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


###############################################################################
# BEGIN RULES


rule all:
	input:
		# expand(
		#	 "sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.tsv",
		#	 simid=simulation_id,
		#	 model_id=MODEL,
		#	 contig=CONTIG,
		#	 rep=REP,
		#	 depth=DEPTH,
		#	 error_rate=ERROR_RATE,
		#	 gc_method=GC_METHODS,
		#	 betavar=BETA_VARS,
		#	 qsbeta=["0_0", "2_5", "2_6"],
		#	 glSpecs=["precise1_gl2", "gl1", "gl2"],
		#	 dogq=[5, 6],
		# ),
		# expand(
		#	 "sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-snp-qs056_intersect.tsv",
		#	 simid=simulation_id,
		#	 model_id=MODEL,
		#	 contig=CONTIG,
		#	 rep=REP,
		#	 depth=DEPTH,
		#	 error_rate=ERROR_RATE,
		#	 gc_method=GC_METHODS,
		#	 betavar=BETA_VARS,
		#	 qsbeta=["0_0", "2_5", "2_6"],
		#	 glSpecs=["precise1_gl2", "gl1", "gl2"],
		#	 dogq=[5, 6],
		# ),
		# expand(
			# "sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.sameSites.bcf",
			# simid=simulation_id,
			# model_id=MODEL,
			# contig=CONTIG,
			# rep=REP,
			# depth=DEPTH,
			# error_rate=ERROR_RATE,
			# gc_method=GC_METHODS,
			# betavar=BETA_VARS,
			# qsbeta=["0_0", "2_5", "2_6"],
			# glSpecs=["precise1_gl2", "gl1", "gl2"],
		# ),
		expand(
			"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ{dogq}.sameSites.tsv",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			gc_method=GC_METHODS,
			betavar=BETA_VARS,
			qsbeta=["0_0", "2_5", "2_6"],
			glSpecs=["precise1_gl2", "gl1", "gl2"],
			dogq=[5, 6],
		),


rule get_pos_tsv:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.tsv",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.tsv",
	shell:
		"""
		{BCFTOOLS} query -f "%CHROM\t%POS\n" {input} > {output}
		"""


rule intersect_pos_tsvs:
	input:
		qs0="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs0_0-snp.tsv",
		qs2_5="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs2_5-snp.tsv",
		qs2_6="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs2_6-snp.tsv",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-snp-qs056_intersect.tsv",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-snp-qs056_intersect.tsv",
	shell:
		"""
		grep -xf <(grep -xf {input.qs0} {input.qs2_5}) {input.qs2_6} > {output}
		"""


rule get_sameSites_bcf:
	input:
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
		tsv="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-snp-qs056_intersect.tsv",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.sameSites.bcf",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.sameSites.bcf",
	shell:
		"""
		({BCFTOOLS} view -T {input.tsv} {input.callgt} -O b -o {output}) 2> {log}
		"""


rule get_genotype_discordance_doGQ_sameSites:
	input:
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.truth.bcf",
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.sameSites.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ{dogq}.sameSites.tsv",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ{dogq}.sameSites.tsv",
	shell:
		"""
		({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output} -doGQ {wildcards.dogq} ) 2> {log}
		"""


