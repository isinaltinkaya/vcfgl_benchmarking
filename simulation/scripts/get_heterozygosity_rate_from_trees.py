#!/opt/software/mamba/23.3.1/bin/python3.10

import tskit


def wattersons_theta(segregating_sites, sample_size):
    """
    Watterson's theta estimator
    """
    a = sum(1.0 / i for i in range(1, int(sample_size)))
    return segregating_sites / a


ts = tskit.load(
    "/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/simulation/sim/sim_vcfgl_2401/model_OutOfAfrica_3G09/contig_chr22/trees/sim_vcfgl_2401-OutOfAfrica_3G09-chr22.trees"
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

#  contig_mutation_rate = 1.29e-08
#
#  print(f"Using contig mutation rate: {contig_mutation_rate}")
#  print(f"Theta_YRI_2 = {theta(Ne_YRI, contig_mutation_rate):.6f}")
#  print(f"Theta_CEU_2 = {theta(Ne_CEU, contig_mutation_rate):.6f}")
#  print(f"Theta_total_2 = {theta(Ne_total, contig_mutation_rate):.6f}")
#


def calculate_heterozygosity(ts):
    return ts.diversity()
    #  heterozygosity = 0
    #  num_sites = 0
#
    #  for site in ts.sites():
        #  num_sites += 1
        #  alleles = [variant.alleles for variant in ts.variants(site=site.id)]
        #  allele_freqs = [sum(variant.genotypes) / len(variant.genotypes) for variant in ts.variants(site=site.id)]
        #  site_heterozygosity = 1 - sum([freq ** 2 for freq in allele_freqs])
        #  heterozygosity += site_heterozygosity
#
    #  return heterozygosity / num_sites if num_sites > 0 else 0
#




print(f"{calculate_heterozygosity(ts):.6f}")



