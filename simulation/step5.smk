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

QSBETAS = ["0_0", "2_5", "2_6", "2_7"]
VCFGLPARS = ["gl1", "gl2", "precise1_gl2", "platform1_gl1", "platform1_gl2"]

BCFTOOLS_MCALL_PRIOR = list(config["bcftools_mcall_prior"].keys())

GC_METHODS = [
    f"{mode}{mcallPriorID}"
    for mode in ["mcall_all_prior", "mcall_perpop_prior"]
    for mcallPriorID in BCFTOOLS_MCALL_PRIOR
] + ["hardcall"]

###############################################################################
# BEGIN RULES


rule all:
    input:
        expand(
            "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-sumSamples_meanReps_all.RData",
            simid=simulation_id,
            model_id=MODEL,
            gc_method=GC_METHODS,
        ),
        # expand(
        #     "sim/{simid}/model_{model_id}/contig_{contig}/gt_discordance/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.tidy.tsv",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     contig=CONTIG,
        #     rep=REP,
        #     depth=DEPTH,
        #     error_rate=ERROR_RATE,
        #     vcfgl_qsbeta=["00", "25", "26", "27"],
        #     vcfgl_gl=["1", "2", "2p"],
        #     vcfgl_platform=[0, 1],
        #     gc_method=GC_METHODS,
        # ),
        # expand(
        #     "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-{gc_method}.tsv",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     contig=CONTIG,
        #     rep=REP,
        #     depth=DEPTH,
        #     error_rate=ERROR_RATE,
        #     vcfgl_qsbeta=["00", "25", "26", "27"],
        #     vcfgl_gl=["1", "2", "2p"],
        #     vcfgl_platform=[0, 1],
        #     gc_method=GC_METHODS,
        # ),
        # expand(
        #     "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-compare_gl12_platform0_qsbeta0567.RData",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     gc_method=GC_METHODS,
        # ),
        # expand(
        #     "sim/{simid}/model_{model_id}/stats/nSNPs/{simid}-{model_id}-{gc_method}-nSNPs.tsv",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     gc_method=GC_METHODS,
        # ),
        # expand(
        #     "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-compare_gl12_platform0_qsbeta0567-R_sumSamplesReps.RData",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     gc_method=GC_METHODS,
        # ),
        # expand(
        #     "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-compare_gl12_platform0_qsbeta0567-R_full.RData",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     gc_method=GC_METHODS,
        # ),
        # expand(
        #     "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-compare_gl12_platform0_qsbeta0567-R_sumSamplesMeanReps.RData",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     gc_method=GC_METHODS,
        # ),


# rule tidy_get_genotype_discordance_doGQ:
#     input:
#         "sim/{simid}/model_{model_id}/contig_{contig}/gt_discordance/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.tsv",
#     output:
#         "sim/{simid}/model_{model_id}/contig_{contig}/gt_discordance/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.tidy.tsv",
#     shell:
#         """
#         awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v VCFGLGL={wildcards.vcfgl_gl} -v VCFGLQSBETA={wildcards.vcfgl_qsbeta} -v VCFGLPLATFORM={wildcards.vcfgl_platform} -v GCMETHOD={wildcards.gc_method} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,VCFGLGL,VCFGLQSBETA,VCFGLPLATFORM,GCMETHOD}}' {input} > {output}
#         """


# rule collect_tidy_get_genotype_discordance_doGQ:
#     input:
#         expand(
#             "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gt_discordance/gcMethod_{{gc_method}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{{gc_method}}.tidy.tsv",
#             contig=CONTIG,
#             rep=REP,
#             depth=DEPTH,
#             error_rate=ERROR_RATE,
#             vcfgl_qsbeta=["00", "25", "26", "27"],
#             vcfgl_gl=["1", "2", "2p"],
#             vcfgl_platform=[0, 1],
#         ),
#     output:
#         "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-{gc_method}.tsv",
#     shell:
#         """
#         cat {input} > {output}
#         """


# rule collect_tidy_nsnps_stats:
#     input:
#         expand(
#             "sim/{{simid}}/model_{{model_id}}/contig_{contig}/stats/nSNPs/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{{gc_method}}.tidy.tsv",
#             contig=CONTIG,
#             rep=REP,
#             depth=DEPTH,
#             error_rate=ERROR_RATE,
#             vcfgl_qsbeta=["00", "25", "26", "27"],
#             vcfgl_gl=["1", "2", "2p"],
#             vcfgl_platform=[0, 1],
#         ),
#     output:
#         "sim/{simid}/model_{model_id}/stats/nSNPs/{simid}-{model_id}-{gc_method}-nSNPs.tsv",
#     shell:
#         """
#         cat {input} > {output}
#         """


# rule save_doGQ_as_RData_sumSamplesReps:
#     input:
#         nSNPs="sim/{simid}/model_{model_id}/stats/nSNPs/{simid}-{model_id}-{gc_method}-nSNPs.tsv",
#         tsvs=expand(
#             "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gt_discordance/gcMethod_{{gc_method}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{{gc_method}}.tidy.tsv",
#             contig=CONTIG,
#             rep=REP,
#             depth=DEPTH,
#             error_rate=ERROR_RATE,
#             vcfgl_qsbeta=["00", "25", "26", "27"],
#             vcfgl_gl=["1", "2"],
#             vcfgl_platform=[0],
#         ),
#     output:
#         "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-compare_gl12_platform0_qsbeta0567-R_sumSamplesReps.RData",
#     shell:
#         """
#         Rscript scripts/save_doGQ_as_RData_sumSamplesReps.R {input.nSNPs} {output} {input.tsvs}
#         """


# rule save_doGQ_as_RData_sumSamplesMeanReps:
#     input:
#         nSNPs="sim/{simid}/model_{model_id}/stats/nSNPs/{simid}-{model_id}-{gc_method}-nSNPs.tsv",
#         tsvs=expand(
#             "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gt_discordance/gcMethod_{{gc_method}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{{gc_method}}.tidy.tsv",
#             contig=CONTIG,
#             rep=REP,
#             depth=DEPTH,
#             error_rate=ERROR_RATE,
#             vcfgl_qsbeta=["00", "25", "26", "27"],
#             vcfgl_gl=["1", "2"],
#             vcfgl_platform=[0],
#         ),
#     output:
#         "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-compare_gl12_platform0_qsbeta0567-R_sumSamplesMeanReps.RData",
#     shell:
#         """
#         Rscript scripts/save_doGQ_as_RData_sumSamplesMeanReps.R {input.nSNPs} {output} {input.tsvs}
#         """


# rule save_doGQ_as_RData:
#     input:
#         nSNPs="sim/{simid}/model_{model_id}/stats/nSNPs/{simid}-{model_id}-{gc_method}-nSNPs.tsv",
#         tsvs=expand(
#             "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gt_discordance/gcMethod_{{gc_method}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{{gc_method}}.tidy.tsv",
#             contig=CONTIG,
#             rep=REP,
#             depth=DEPTH,
#             error_rate=ERROR_RATE,
#             vcfgl_qsbeta=["00", "25", "26", "27"],
#             vcfgl_gl=["1", "2", "2p"],
#             vcfgl_platform=[0,1],
#         ),
#     output:
#         "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-compare_gl12_platform0_qsbeta0567-R_full.RData",
#     shell:
#         """
#         Rscript scripts/save_doGQ_as_RData.R {input.nSNPs} {output} {input.tsvs}
#         """


# rule save_doGQ_as_RData:
#     input:
#         expand(
#             "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gt_discordance/gcMethod_{{gc_method}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{{gc_method}}.tidy.tsv",
#             contig=CONTIG,
#             rep=REP,
#             depth=DEPTH,
#             error_rate=ERROR_RATE,
#             vcfgl_qsbeta=["00", "25", "26", "27"],
#             vcfgl_gl=["1", "2", "2p"],
#             vcfgl_platform=[0, 1],
#         ),
#     output:
#         "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-compare_all.RData",
#     shell:
#         """
#         Rscript scripts/save_doGQ_as_RData.R {output} {input}
#         """


rule save_doGQ_as_RData_sumSamplesMeanReps:
    input:
        expand(
            "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gt_discordance/gcMethod_{{gc_method}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{{gc_method}}.tidy.tsv",
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            vcfgl_qsbeta=["00", "25", "26", "27"],
            vcfgl_gl=["1", "2", "2p"],
            vcfgl_platform=[0, 1],
        ),
    output:
        "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance/{simid}-{model_id}-gcMethod{gc_method}-sumSamples_meanReps_all.RData",
    shell:
        """
        Rscript scripts/save_doGQ_as_RData_sumSamplesMeanReps.R {output} {input}
        """
