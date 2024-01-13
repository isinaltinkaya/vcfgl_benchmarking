#!/opt/software/mamba/23.3.1/bin/python3.10

import tskit


def wattersons_theta(segregating_sites, sample_size):
    """
    Watterson's theta estimator
    """
    a = sum(1.0 / i for i in range(1, int(sample_size)))
    return segregating_sites / a


ts = tskit.load(
    "/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/simulation/sim/sim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/trees/sim_vcfgl_2312-OutOfAfrica_3G09-chr22.trees"
)

# We simulate 50 inds from YRI and 50 inds from CEU
simulated_populations = ["YRI", "CEU"]
npops = len(simulated_populations)

# get the sample set from each pop
sample_sets = [ts.samples(population=pop_id) for pop_id in range(npops)]
# print(sample_sets)

seg_sites_per_pop = {
    pop_id: ts.segregating_sites(ts.samples(pop_id), span_normalise=True)
    for pop_id in range(npops)
}
# print(seg_sites_per_pop)

theta_per_pop = {}
for pop_id, seg_sites in seg_sites_per_pop.items():
    sample_size = sum(1 for i in ts.samples() if ts.node(i).population == pop_id)
    theta_per_pop[pop_id] = wattersons_theta(seg_sites, sample_size)

popnames = ["YRI", "CEU"]

for pop_id, theta in theta_per_pop.items():
    print(f"\hat{{Theta_{popnames[pop_id]}}} = {theta:.6f}")


theta_total = wattersons_theta(ts.segregating_sites(), ts.num_samples)
print(f"\hat{{Theta_total}} = {theta_total:.6f}")


mutation_rate = 2.35e-8
Ne_YRI = 12300
# CEU pop. size after EU/AS divergence
Ne_CEU = 1000
# CHB pop. size after EU/AS divergence
Ne_CHB = 510

Ne_total = Ne_YRI + Ne_CEU + Ne_CEU


def theta(Ne, mu):
    return 4 * Ne * mu


print(f"Using model mutation rate: {mutation_rate}")
print(f"Theta_YRI = {theta(Ne_YRI, mutation_rate):.6f}")
print(f"Theta_CEU = {theta(Ne_CEU, mutation_rate):.6f}")
print(f"Theta_total = {theta(Ne_total, mutation_rate):.6f}")

contig_mutation_rate = 1.29e-08

print(f"Using contig mutation rate: {contig_mutation_rate}")
print(f"Theta_YRI_2 = {theta(Ne_YRI, contig_mutation_rate):.6f}")
print(f"Theta_CEU_2 = {theta(Ne_CEU, contig_mutation_rate):.6f}")
print(f"Theta_total_2 = {theta(Ne_total, contig_mutation_rate):.6f}")


# # print("\n\n\n\n\n")
# # outf = open("delme.vcf", "w")
# # ts.write_vcf(outf, ploidy=2, individual_names=indv_names)
# # outf.close()

# import msprime
# import stdpopsim
# import numpy as np
# import sys
# import yaml

# ts = msprime.sim_ancestry(samples=3, ploidy=2, sequence_length=10, random_seed=2)
# ts = msprime.sim_mutations(ts, rate=0.1, random_seed=2)
# ts.write_vcf(sys.stdout, individual_names=["A", "B", "C"])


# for var in ts.variants():
#     #print if variant is not biallelic
#     if len(var.alleles) > 2:
#         print(var.site.position, var.alleles, var.genotypes, sep="\t")


# print(ts)

# configfile="/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/simulation/config/sim_vcfgl_2401.yaml"
# configfile="config/sim_vcfgl_2401.yaml"
# with open(
#     configfile,
#     "r",
# ) as ymlfile:
#     config = yaml.safe_load(ymlfile)

# simulation_id = config["simulation_id"]
# model_id = config["model"]
# contig_id = config["contig"]


# DEF_POPS = config["def_pops"]
# ploidy = 2

# haplo_list = []
# indv_names = []


# DEF_POPS = {"YRI": 50, "CEU": 50, "CHB": 0}
# ploidy = 2
# pops = ["YRI", "CEU"]

# haplo_list = []
# indv_names = []

# for i, (key, value) in enumerate(DEF_POPS.items()):
#     haplo_list.append(value * ploidy)
#     for ind in range(value):
#         # Using PLINK-like format: <Family-ID>_<Individual-ID>
#         # to store <Population-ID>_<Individual-ID>
#         indv_names.append(f"pop{key}_ind{str(ind+1)}")

# # print(haplo_list)
# print(indv_names)

# samples = DEF_POPS

# species = stdpopsim.get_species("HomSap")
# model = species.get_demographic_model(model_id)
# contig = species.get_contig(contig_id)
# # contig = species.get_contig("chr22", mutation_rate=model.mutation_rate)

# # samples = model.get_samples(*haplo_list)
# engine = stdpopsim.get_engine("msprime")
# seedval = np.random.randint(2**20) + 42

# samples = model.get_samples(*haplo_list)

# seedval=829581

# ts = engine.simulate(model, contig, samples, seed=seedval)
# ts2 = engine.simulate(model, contig, samples, seed=seedval, msprime_model="dtwf")

# ts3=tskit.load("sim/prevsim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/trees/sim_vcfgl_2312-OutOfAfrica_3G09-chr22.trees")


# print(ts3)

# ts.write_vcf(sys.stdout)

# with open("delme.ts", "w") as tsout:
#     ts.dump(tsout)

# with open("ts1.ts", "w") as tsout:
#     ts.dump(tsout)

# with open("ts2.ts", "w") as tsout:
#     ts2.dump(tsout)
# with open("ts3.ts", "w") as tsout:
#     ts3.dump(tsout)
# with open("ts1.vcf", "w") as vcfout:
#     ts.write_vcf(vcfout)
# with open("ts2.vcf", "w") as vcfout:
#     ts2.write_vcf(vcfout)
# with open("ts3.vcf", "w") as vcfout:
#     ts3.write_vcf(vcfout)
# TODO!!
# NOTE
# /home/pfs488/.local/lib/python3.10/site-packages/stdpopsim/engines.py:120: UserWarning: The demographic model has mutation rate 2.35e-08, but this simulation used the contig's mutation rate 1.29e-08. Diversity levels may be different than expected for this species.
