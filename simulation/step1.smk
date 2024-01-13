###############################################################################
# Step 1: Simulate data for vcfgl benchmarking
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


simulation_id = config["simulation_id"]


MODEL = config["model"]
CONTIG = config["contig"]

species = stdpopsim.get_species("HomSap")

DEF_POPS = config["def_pops"]
ploidy = 2

haplo_list = []
indv_names = []

for i, (key, value) in enumerate(DEF_POPS.items()):
    haplo_list.append(value * ploidy)
    for ind in range(value):
        # Using PLINK-like format: <Family-ID>_<Individual-ID>
        # to store <Population-ID>_<Individual-ID>
        indv_names.append(f"pop{key}_ind{str(ind+1)}")

###############################################################################
# BEGIN RULES


rule all:
    input:
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
        ),
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
        ),


# ________________________________________________________________________________
# BEGIN SIMULATION


rule simulation:
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
    params:
        logfile="sim/{simid}/logs/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
    run:
        model = species.get_demographic_model(wildcards.model_id)
        contig = species.get_contig(wildcards.contig)
        samples = model.get_samples(*haplo_list)
        engine = stdpopsim.get_engine("msprime")
        seedval = np.random.randint(2**20) + 42
        with open(params.logfile, "w") as logout:
            print(
                "Contig index: "
                + str(int(str(wildcards.contig)[3:5]))
                + ", chr: "
                + str(wildcards.contig),
                file=logout,
            )
            print("Seed value: " + str(seedval), file=logout)
        ts = engine.simulate(model, contig, samples, seed=seedval, msprime_model="dtwf")
        with open(output[0], "w") as tsout:
            ts.dump(tsout)


## Convert tree sequence to VCF file
# Using legacy format to avoid multiple instances of sites
rule tree_to_vcf:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf",
    run:
        ts = tskit.load(input[0])
        with open(output[0], "w") as vcfout:
            ts.write_vcf(
                output=vcfout,
                ploidy=ploidy,
                contig_id=str(wildcards.contig),
                position_transform="legacy",
                individual_names=indv_names,
            )


# TODO
rule print_pop_inds_list:
    output:
        config["pop_inds_list"],
    run:
        with open(output[0], "w") as pop_inds_list:
            for i, indv in enumerate(indv_names):
                print(indv + "\tpop" + str(i // ploidy + 1), file=pop_inds_list)
