###############################################################################
# Step 5: Collect and plot genotype discordance metrics
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

MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]

DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]

GC_METHODS = ["genotype_calling", "genotype_calling_perpop"]

QSBETAS = ["0_0", "2_5", "2_6", "2_7"]
GLSPECS = ["gl1", "gl2", "precise1_gl2"]


###############################################################################
# BEGIN RULES


rule all:
    input:
        expand(
			"sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}.tsv",
			# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.tsv.tidy",
            qsbeta=["0_0", "2_5", "2_6", "2_7"],
            # "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}.tsv",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
			gc_method=GC_METHODS,
            glSpecs=GLSPECS,
        ),
        # expand(
			# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.tsv.tidy",
            # qsbeta=["0_0", "2_5", "2_6", "2_7"],
            # simid=simulation_id,
            # model_id=MODEL,
            # contig=CONTIG,
            # rep=REP,
            # depth=DEPTH,
            # error_rate=ERROR_RATE,
			# gc_method=GC_METHODS,
            # glSpecs=GLSPECS,
        # ),
        # expand(
            # "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-hard_genotype_calling.tsv",
            # simid=simulation_id,
            # model_id=MODEL,
            # contig=CONTIG,
            # rep=REP,
            # depth=DEPTH,
            # error_rate=ERROR_RATE,
            # glSpecs=GLSPECS,
        # ),
        # expand(
        #     "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}-doGQ7.tsv",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     contig=CONTIG,
        #     rep=REP,
        #     depth=DEPTH,
        #     error_rate=ERROR_RATE,
        #     gc_method=GC_METHODS,
        #     glSpecs=GLSPECS,
        # ),
        # expand(
        # "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}-doGQ{dogq}.tsv",
        # simid=simulation_id,
        # model_id=MODEL,
        # contig=CONTIG,
        # rep=REP,
        # depth=DEPTH,
        # error_rate=ERROR_RATE,
        # gc_method=GC_METHODS,
        # dogq=[6],
        # glSpecs=["precise1_gl2", "gl1", "gl2"],
        # ),



rule tidy_get_genotype_discordance_snp:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}-snp.tsv",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}-snp.tsv.tidy",
	shell:
		"""
		awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v BETAVAL={wildcards.betaval} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,BETAVAL}}' {input} > {output}
		"""


rule collect_get_genotype_discordance_qs_snp:
	input:
		expand(
			"sim/{{simid}}/model_{{model_id}}/contig_{contig}/gc_evaluation/genotype_discordance/{{gc_method}}_{{glSpecs}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.tsv.tidy",
			contig=CONTIG,
			depth=DEPTH,
			rep=REP,
			error_rate=ERROR_RATE,
			qsbeta=["0_0", "2_5", "2_6", "2_7"],
		),
	output:
		"sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}-snp.tsv"
	shell:
		"""
		cat {input} > {output}
		"""


rule tidy_get_genotype_discordance:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}.tsv",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}.tsv.tidy",
	shell:
		"""
		awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v BETAVAL={wildcards.betaval} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,BETAVAL}}' {input} > {output}
		"""


rule collect_get_genotype_discordance_qs:
	input:
		expand(
			"sim/{{simid}}/model_{{model_id}}/contig_{contig}/gc_evaluation/genotype_discordance/{{gc_method}}_{{glSpecs}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv.tidy",
			contig=CONTIG,
			depth=DEPTH,
			rep=REP,
			error_rate=ERROR_RATE,
			qsbeta=["0_0", "2_5", "2_6", "2_7"],
		),
	output:
		"sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}.tsv"
	shell:
		"""
		cat {input} > {output}
		"""


# rule tidy_get_genotype_discordance_hardcall:
    # input:
        # "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/hard_genotype_calling_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}.tsv",
    # output:
        # "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/hard_genotype_calling_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}.tsv.tidy",
    # shell:
        # """
        # awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v BETAVAL={wildcards.betaval} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,BETAVAL}}' {input} > {output}
        # """
#
#
# rule collect_get_genotype_discordance_qs_hardcall:
    # input:
        # expand(
            # "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gc_evaluation/genotype_discordance/hard_genotype_calling_{{glSpecs}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv.tidy",
            # contig=CONTIG,
            # depth=DEPTH,
            # rep=REP,
            # error_rate=ERROR_RATE,
            # qsbeta=["0_0", "2_5", "2_6", "2_7"],
        # ),
    # output:
        # "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-hard_genotype_calling.tsv",
    # shell:
        # """
        # cat {input} > {output}
        # """
#
#
# rule tidy_get_genotype_discordance_doGQ:
    # input:
        # "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}-snp_doGQ{dogq}.tsv",
    # output:
        # "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}-snp_doGQ{dogq}.tsv.tidy",
    # shell:
        # """
        # awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v BETAVAL={wildcards.betaval} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,BETAVAL}}' {input} > {output}
        # """
#
#
# rule tidy_get_genotype_discordance_doGQ7:
    # input:
        # "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}_doGQ7.tsv",
    # output:
        # "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betaval}_doGQ7.tsv.tidy",
    # shell:
        # """
        # awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v BETAVAL={wildcards.betaval} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,BETAVAL}}' {input} > {output}
        # """
#
#
# rule collect_get_genotype_discordance_qs_doGQ:
    # input:
        # expand(
            # "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gc_evaluation/genotype_discordance/{{gc_method}}_{{glSpecs}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ{{dogq}}.tsv.tidy",
            # contig=CONTIG,
            # depth=DEPTH,
            # rep=REP,
            # error_rate=ERROR_RATE,
            # qsbeta=["0_0", "2_5", "2_6"],
        # ),
    # output:
        # "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}-doGQ{dogq}.tsv",
    # shell:
        # """
        # cat {input} > {output}
        # """
#
#
# rule collect_get_genotype_discordance_qs_doGQ7:
    # input:
        # expand(
            # "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gc_evaluation/genotype_discordance/{{gc_method}}_{{glSpecs}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}_doGQ7.tsv.tidy",
            # contig=CONTIG,
            # depth=DEPTH,
            # rep=REP,
            # error_rate=ERROR_RATE,
            # qsbeta=["0_0", "2_5", "2_6", "2_7"],
        # ),
    # output:
        # "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_{glSpecs}/{simid}-{model_id}-{gc_method}-doGQ7.tsv",
    # shell:
        # """
        # cat {input} > {output}
        # """
