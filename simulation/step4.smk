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


GET_GT_DISCORDANCE = config["tools"]["gtDiscordance"]


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
        expand("sim/{simid}/model_{model_id}/contig_{contig}/stats/nSNPs.txt",simid=simulation_id,model_id=MODEL,contig=CONTIG)
        # expand(
        #     "sim/{simid}/model_{model_id}/contig_{contig}/gt_discordance/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.tsv",
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
        #     "sim/{simid}/model_{model_id}/stats/nSites_non0dp/{simid}-{model_id}_nSitesNon0DP.csv",
        #     simid=simulation_id,
        #     model_id=MODEL,
        # ),
        # expand(
            # "sim/{simid}/model_{model_id}/contig_{contig}/stats/bcfstats/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.txt",
            # simid=simulation_id,
            # model_id=MODEL,
            # contig=CONTIG,
            # rep=REP,
            # depth=DEPTH,
            # error_rate=ERROR_RATE,
            # vcfgl_qsbeta=["00", "25", "26", "27"],
            # vcfgl_gl=["1", "2", "2p"],
            # vcfgl_platform=[0, 1],
            # gc_method=GC_METHODS,
        # ),
        # expand(
            # "sim/{simid}/model_{model_id}/contig_{contig}/stats/nSNPs/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.tidy.tsv",
            # simid=simulation_id,
            # model_id=MODEL,
            # contig=CONTIG,
            # rep=REP,
            # depth=DEPTH,
            # error_rate=ERROR_RATE,
            # vcfgl_qsbeta=["00", "25", "26", "27"],
            # vcfgl_gl=["1", "2", "2p"],
            # vcfgl_platform=[0, 1],
            # gc_method=GC_METHODS,
        # ),
        # expand(
            # "sim/{simid}/model_{model_id}/contig_{contig}/gt_discordance/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.snp.tsv",
            # simid=simulation_id,
            # model_id=MODEL,
            # contig=CONTIG,
            # rep=REP,
            # depth=DEPTH,
            # error_rate=ERROR_RATE,
            # vcfgl_qsbeta=["00", "25", "26", "27"],
            # vcfgl_gl=["1", "2", "2p"],
            # vcfgl_platform=[0, 1],
            # gc_method=GC_METHODS,
        # ),


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


# rule get_perSample_non0depth_stats:
#     input:
#         "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform0_gl2_qs00/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform0_gl2_qs00.bcf",
#     output:
#         "sim/{simid}/model_{model_id}/contig_{contig}/stats/perSample_nSitesNon0DP/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}.ind{indv}.csv",
#     shell:
#         """
#         # sample id, rep, depth, nSitesNon0DP
#         printf "{wildcards.indv},{wildcards.rep},{wildcards.depth},$({BCFTOOLS} view -s {wildcards.indv} {input} | {BCFTOOLS} filter -i 'FMT/DP>0' | {BCFTOOLS} view -H | wc -l)\n" > {output}
#         """


# rule merge_perSample_non0depth_stats:
#     input:
#         expand(
#             "sim/{{simid}}/model_{{model_id}}/contig_{contig}/stats/perSample_nSitesNon0DP/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}.ind{indv}.csv",
#             contig=CONTIG,
#             depth=DEPTH,
#             rep=REP,
#             error_rate=ERROR_RATE,
#             indv=indv_names,
#         ),
#     output:
#         "sim/{simid}/model_{model_id}/stats/nSites_non0dp/{simid}-{model_id}_nSitesNon0DP.csv",
#     shell:
#         """
#         (echo "Sample,Rep,Depth,nSitesNon0DP";
#         for i in {input}; do
#         cat $i;
#         done) > {output}
#         """


# rule get_genotype_discordance_hardCall:
# input:
# truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.truth.bcf",
# callgt="sim/{simid}/model_{model_id}/contig_{contig}/hardcall_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
# output:
# "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/hardcall_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
# log:
# "sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/hardcall_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
# shell:
# """
# ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output}) 2> {log}
# """

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


# rule get_genotype_discordance_otherGcMethods:
# 	 input:
# 		 truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qsbeta}.truth.bcf",
# 		 callgt="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
# 	 output:
# 		 "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
# 	 log:
# 		 "sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{glSpecs}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
# 	 shell:
# 		 """
# 		 ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output}) 2> {log}
# 		 """


# rule get_genotype_discordance_gcMethods:
#     input:
#         truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_{vcfglpars}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.truth.bcf",
#         callgt="sim/{simid}/model_{model_id}/contig_{contig}/{gc_method}_{vcfglpars}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.bcf",
#     output:
#         "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{vcfglpars}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
#     log:
#         "sim/{simid}/logs/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_{vcfglpars}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}.tsv",
#     shell:
#         """
#         ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output}) 2> {log}
#         """


def get_dogq_for_gc_method(wildcards):
    if wildcards.gc_method == "hardcall":
        return "8"
    elif wildcards.gc_method.startswith("mcall"):
        return "7"
    else:
        raise ValueError(f"gc_method {wildcards.gc_method} not recognized")


rule get_genotype_discordance_forSNPs_doGQ:
    input:
        truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.truth.bcf",
        callgt="sim/{simid}/model_{model_id}/contig_{contig}/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/gt_discordance/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.snp.tsv",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/gt_discordance/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.snp.tsv",
    params:
        dogqval=get_dogq_for_gc_method,
    shell:
        """
        {BCFTOOLS} view -v snps {input.callgt} -O b -o {input.callgt}_snp.bcf;
        ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt}_snp.bcf -o {output} -doGQ {params.dogqval} ) 2> {log}
        """


rule get_genotype_discordance_doGQ:
    input:
        truthgt="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.truth.bcf",
        callgt="sim/{simid}/model_{model_id}/contig_{contig}/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/gt_discordance/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.tsv",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/gt_discordance/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.tsv",
    params:
        dogqval=get_dogq_for_gc_method,
    shell:
        """
        ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output} -doGQ {params.dogqval} ) 2> {log}
        """


rule get_bcfstats:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/gcMethod_{gc_method}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/stats/bcfstats/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.txt",
    shell:
        """
        {BCFTOOLS} stats {input} > {output}
        """


rule get_n_snps_and_tidy:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/stats/bcfstats/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.txt",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/stats/nSNPs/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.tidy.tsv",
    shell:
        """
        awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v VCFGLGL={wildcards.vcfgl_gl} -v VCFGLQSBETA={wildcards.vcfgl_qsbeta} -v VCFGLPLATFORM={wildcards.vcfgl_platform} -v GCMETHOD={wildcards.gc_method} 'BEGIN{{FS="\\t";OFS="\\t"}} /^SN/ && /number of SNPs/ {{split($0, a, ":"); split(a[2], b, "\\t"); print b[2],REP,DEPTH,VCFGLGL,VCFGLQSBETA,VCFGLPLATFORM,GCMETHOD}}' {input} > {output}
        """


rule collect_get_n_snps_and_tidy:
    input:
        expand("sim/{simid}/model_{model_id}/contig_{contig}/stats/nSNPs/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_{gc_method}.tidy.tsv",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            vcfgl_qsbeta=["00", "25", "26", "27"],
            vcfgl_gl=["1", "2"],
            vcfgl_platform=[0, 1],
            gc_method=GC_METHODS)
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/stats/nSNPs.txt"
    shell:
        """
		cat {input} > {output}
        """



