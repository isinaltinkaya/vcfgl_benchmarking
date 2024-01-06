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


simulation_id = config["simulation_id"]



MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]
DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]


POPS = config["pops"]

GLMODELS = [1, 2]

ANGSD = config["tools"]["angsd"]
BCFTOOLS = config["tools"]["bcftools"]


###############################################################################
# BEGIN RULES


rule all:
    input:
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_perpop_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
            glSpecs=["precise1_gl2", "gl1", "gl2"],
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            qsbeta=["0_0", "2_5", "2_6"],
            population=POPS,
        ),
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
            glSpecs=["precise1_gl2", "gl1", "gl2"],
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            qsbeta=["0_0", "2_5", "2_6"],
            population=POPS,
        ),


# ###############################################################################


rule bcftools_call_snp_and_genotypes:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
    params:
        prior=config["bcftools_call_m_prior"],
    shell:
        """
         ({BCFTOOLS} call --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 -mv --prior {params.prior} -O b -o {output} {input} )2> {log}
         """


###############################################################################


rule bcftools_call_perpop_G:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_perpop_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
    params:
        prior=config["bcftools_call_m_prior"],
        popinds_list="sim/sim_vcfgl_2312/sim_vcfgl_2312_samples_pops.tsv",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling_perpop_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp.bcf",
    shell:
        """
        ({BCFTOOLS} call -G {params.popinds_list} --annotate FORMAT/GQ,INFO/PV4 -mv --prior {params.prior} -O b -o {output} {input} )2> {log}
        """
