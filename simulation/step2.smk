###############################################################################
# Step 2: Run vcfgl simulation with qs error
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

simulation_id = "sim_vcfgl_2312"


configfile: "config/" + simulation_id + ".yaml"


MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]
DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]

BETA_VARS = config["beta_variance_values_neg_e"]
qsbeta = ["2_" + str(x) for x in BETA_VARS]

VCFGL = config["tools"]["vcfgl"]

###############################################################################
# BEGIN RULES


GLMODELS = [1, 2]


rule all:
    input:
        # expand(
            # "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.bcf",
            # simid=simulation_id,
            # model_id=MODEL,
            # contig=CONTIG,
            # rep=REP,
            # depth=DEPTH,
            # error_rate=ERROR_RATE,
            # qsbeta=["0_0", "2_5", "2_6", "2_7"],
            # glModel=1,
        # ),
        expand("sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_precise1_gl2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.bcf",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            qsbeta=["0_0", "2_5", "2_6", "2_7"],
        ),


def get_error_params(wildcards):
	if wildcards.qsbeta == "0_0":
		return " --error-qs 0 "
	else:
		return str(
			" --error-qs "
			+ wildcards.qsbeta.split("_")[0]
			+ " --beta-variance 1e-"
			+ wildcards.qsbeta.split("_")[1]
		)


# rule simulate_vcfgl_var_qs_error:
    # input:
        # "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf",
    # output:
        # "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.bcf",
    # params:
        # prefix="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}",
        # random_seed=lambda wildcards: str(int(wildcards.rep) + 1),
        # errparam=get_error_params,
    # log:
        # "sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.bcf",
    # shell:
        # """
        # ({VCFGL} -i {input} \
                # -o {params.prefix} \
                # -O b \
                # --depth {wildcards.depth} \
                # --error-rate {wildcards.error_rate}  \
                # {params.errparam} \
                # -addQS 1 \
                # -explode 0 \
                # -addInfoAD 1 \
                # -addFormatAD 1 \
                # -addFormatDP 1 \
                # -addPL 1 \
                # -GL {wildcards.glModel} \
                # -printTruth 1 \
                # -printPileup 1 \
                # --rm-invar-sites 3 \
                # --rm-empty-sites 1 \
                # -doUnobserved 1 \
                # --seed {params.random_seed} ) 2> {log}
        # """
#
#
rule simulate_vcfgl_var_qs_error_precisegl1:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_precise1_gl2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.bcf",
    params:
        prefix="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_precise1_gl2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}",
        random_seed=lambda wildcards: str(int(wildcards.rep) + 1),
        errparam=get_error_params,
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_precise1_gl2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.bcf",
    shell:
        """
        ({VCFGL} -i {input} \
                -o {params.prefix} \
                -O b \
                --depth {wildcards.depth} \
                --error-rate {wildcards.error_rate}  \
                {params.errparam} \
                -addQS 1 \
                -explode 0 \
                -addInfoAD 1 \
                -addFormatAD 1 \
                -addFormatDP 1 \
                -addPL 1 \
                -GL 2 \
                -printTruth 1 \
                -printPileup 1 \
                --rm-invar-sites 3 \
                --rm-empty-sites 1 \
				--precise-gl 1 \
                -doUnobserved 1 \
                --seed {params.random_seed} ) 2> {log}
        """
