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


QSBETAS = ["0_0", "2_5", "2_6", "2_7"]
GLSPECS = ["gl1", "gl2", "precise1_gl2"]


###############################################################################
# BEGIN RULES


rule all:
	input:
		expand(
			"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			gc_method=GC_METHODS,
			betavar=BETA_VARS,
			qsbeta=QSBETAS,
			glSpecs=GLSPECS,
		),
		# expand(
			# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.tsv",
			# simid=simulation_id,
			# model_id=MODEL,
			# contig=CONTIG,
			# rep=REP,
			# depth=DEPTH,
			# error_rate=ERROR_RATE,
			# gc_method=GC_METHODS,
			# betavar=BETA_VARS,
			# qsbeta=QSBETAS,
			# glSpecs=GLSPECS,
		# ),
		# expand(
			# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/hard_genotype_calling_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
			# simid=simulation_id,
			# model_id=MODEL,
			# contig=CONTIG,
			# rep=REP,
			# depth=DEPTH,
			# error_rate=ERROR_RATE,
			# gc_method="hard_genotype_call",
			# betavar=BETA_VARS,
			# qsbeta=QSBETAS,
			# glSpecs=GLSPECS,
		# ),
		# expand(
		#	 "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}_doGQ7.tsv",
		#	 simid=simulation_id,
		#	 model_id=MODEL,
		#	 contig=CONTIG,
		#	 rep=REP,
		#	 depth=DEPTH,
		#	 error_rate=ERROR_RATE,
		#	 gc_method=GC_METHODS,
		#	 betavar=BETA_VARS,
		#	 qsbeta=QSBETAS,
		#	 glSpecs=GLSPECS,
		# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ{dogq}.tsv",
		# 		simid=simulation_id,
		# 		model_id=MODEL,
		# 		contig=CONTIG,
		# 		rep=REP,
		# 		depth=DEPTH,
		# 		error_rate=ERROR_RATE,
		# 		gc_method=GC_METHODS,
		# 		betavar=BETA_VARS,
		# 		qsbeta=QSBETAS,
		# 		glSpecs=GLSPECS,
		# 		dogq=[6]),
		#  expand("sim/{simid}/model_{model_id}/stats/nSites_non0dp/{simid}-{model_id}_nSitesNon0DP.csv",
		# 		simid=simulation_id,
		# 		model_id=MODEL)

#
# rule get_genotype_discordance_doGQ:
	# input:
		# truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.truth.bcf",
		# callgt="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ{dogq}.tsv",
	# log:
		# "sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ{dogq}.tsv",
	# shell:
		# """
		# ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output} -doGQ {wildcards.dogq} ) 2> {log}
		# """
#
#
# rule get_genotype_discordance_doGQ7_nonSNP:
	# input:
		# truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.truth.bcf",
		# callgt="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}_doGQ7.tsv",
	# log:
		# "sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}_doGQ7.tsv",
	# shell:
		# """
		# ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output} -doGQ 7 ) 2> {log}
		# """
#
#
# # expand(
# # "sim/{simid}/model_{model_id}/stats/nSites_non0dp/{simid}-{model_id}_nSitesNon0DP.csv",
# # simid=simulation_id,
# # model_id=MODEL,
# # ),
# # N.B. using gl 1; gl model does not matter since with same seed vcfgl poisson samples the same depths
# # validation:
# # $ bcftools view -s popCEU_ind1 sim/sim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/vcfgl_gl1/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d2-e0.002_qs0_0.bcf  | bcftools filter -i 'FMT/DP>0' | bcftools view -H | wc -l
# # 157954
# # $ bcftools view -s popCEU_ind1 sim/sim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/vcfgl_gl2/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d2-e0.002_qs0_0.bcf  | bcftools filter -i 'FMT/DP>0' | bcftools view -H | wc -l
# # 157954
# #
# # N.B. using qsbeta 0_0; qsbeta does not matter since with same seed vcfgl poisson samples the same depths
# # validation:
# # $ bcftools view -s popCEU_ind1 sim/sim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/vcfgl_gl2/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d2-e0.002_qs0_0.bcf  | bcftools filter -i 'FMT/DP>0' | bcftools view -H | wc -l
# # 157954
# # $ bcftools view -s popCEU_ind1 sim/sim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/vcfgl_gl2/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d2-e0.002_qs2_5.bcf  | bcftools filter -i 'FMT/DP>0' | bcftools view -H | wc -l
# # 157954
# rule get_perSample_non0depth_stats:
	# input:
		# "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_gl1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs0_0.bcf",
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/stats/perSample_nSitesNon0DP/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}.ind{indv}.csv",
	# shell:
		# """
		 # # sample id, rep, depth, nSitesNon0DP
		 # printf "{wildcards.indv},{wildcards.rep},{wildcards.depth},$({BCFTOOLS} view -s {wildcards.indv} {input} | {BCFTOOLS} filter -i 'FMT/DP>0' | {BCFTOOLS} view -H | wc -l)\n" > {output}
		 # """
#
#
# rule merge_perSample_non0depth_stats:
	# input:
		# expand(
			# "sim/{{simid}}/model_{{model_id}}/contig_{contig}/stats/perSample_nSitesNon0DP/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}.ind{indv}.csv",
			# contig=CONTIG,
			# depth=DEPTH,
			# rep=REP,
			# error_rate=ERROR_RATE,
			# indv=indv_names,
		# ),
	# output:
		# "sim/{simid}/model_{model_id}/stats/nSites_non0dp/{simid}-{model_id}_nSitesNon0DP.csv",
	# shell:
		# """
		 # (echo "Sample,Rep,Depth,nSitesNon0DP";
		 # for i in {input}; do
			 # cat $i;
		 # done) > {output}
		 # """
#
#
# rule get_genotype_discordance_hardCall:
	# input:
		# truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.truth.bcf",
		# callgt="sim/{simid}/model_{model_id}/contig_{contig}/hard_genotype_calling_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/hard_genotype_calling_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
	# log:
		# "sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/hard_genotype_calling_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
	# shell:
		# """
		# ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output}) 2> {log}
		# """
#
# rule get_genotype_discordance_otherGcMethods_snp:
	# input:
		# truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.truth.bcf",
		# callgt="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.tsv",
	# log:
		# "sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.tsv",
	# shell:
		# """
		# ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output}) 2> {log}
		# """
#
rule get_genotype_discordance_otherGcMethods:
	input:
		truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.truth.bcf",
		callgt="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
	shell:
		"""
		({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output}) 2> {log}
		"""



