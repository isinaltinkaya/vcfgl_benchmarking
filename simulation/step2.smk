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


simulation_id = config["simulation_id"]


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
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.bcf",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            vcfgl_qsbeta=["00", "25", "26", "27"],
            vcfgl_gl=["1", "2", "2p"],
            vcfgl_platform=[0, 1],
        ),


def get_gl_params(wildcards):
    if wildcards.vcfgl_gl == "1":
        return " --precise-gl 0 --gl-model 1 "
    elif wildcards.vcfgl_gl == "2":
        return " --precise-gl 0 --gl-model 2 "
    elif wildcards.vcfgl_gl == "2p":
        return " --precise-gl 1 --gl-model 2 "
    else:
        raise ValueError("Invalid gl model")


def get_error_params(wildcards):
    if wildcards.vcfgl_qsbeta == "00":
        return " --error-qs 0 "
    else:
        return str(
            " --error-qs "
            + wildcards.vcfgl_qsbeta[0]
            + " --beta-variance 1e-"
            + wildcards.vcfgl_qsbeta[1]
        )

#
# rule simulate_vcfgl:
    # input:
        # "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf",
    # output:
        # "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.bcf",
    # params:
        # prefix="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}",
        # random_seed=lambda wildcards: str(int(wildcards.rep) + 1),
        # errparam=get_error_params,
        # glparam=get_gl_params,
    # log:
        # "sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.bcf",
    # shell:
        # """
        # ({VCFGL} -i {input} \
                # -o {params.prefix} \
                # -O b \
                # --depth {wildcards.depth} \
                # --error-rate {wildcards.error_rate}  \
                # {params.errparam} \
                # -addQS 1 \
                # -explode 1 \
                # -addInfoAD 1 \
                # -addFormatAD 1 \
                # -addFormatDP 1 \
                # -addPL 1 \
                # {params.glparam} \
                # -printTruth 1 \
                # -printPileup 1 \
                # --rm-invar-sites 3 \
                # --rm-empty-sites 1 \
                # -doUnobserved 1 \
                # --source 1 \
                # --platform {wildcards.vcfgl_platform} \
                # --seed {params.random_seed} ) 2> {log}
        # """
#

# 240131 update for sim_2402
# --rm-invar-sites 0 \
# --rm-empty-sites 0 \
rule simulate_vcfgl:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.bcf",
    params:
        prefix="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}",
        random_seed=lambda wildcards: str(int(wildcards.rep) + 1),
        errparam=get_error_params,
        glparam=get_gl_params,
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.bcf",
    shell:
        """
        ({VCFGL} -i {input} \
                -o {params.prefix} \
                -O b \
                --depth {wildcards.depth} \
                --error-rate {wildcards.error_rate}  \
                {params.errparam} \
                -addQS 1 \
                -explode 1 \
                -addInfoAD 1 \
                -addFormatAD 1 \
                -addFormatDP 1 \
                -addPL 1 \
                {params.glparam} \
                -printTruth 1 \
                -printPileup 1 \
                --rm-invar-sites 0 \
                --rm-empty-sites 0 \
                -doUnobserved 1 \
                --source 1 \
                --platform {wildcards.vcfgl_platform} \
                --seed {params.random_seed} ) 2> {log}
        """

