#include <stdio.h>
#include <htslib/vcf.h>

/*
 * Macro:[AT]
 * inject the file and line info as string
 */
#define STRINGIFY(x) #x
#define ASSTR(x) STRINGIFY(x)
#define AT __FILE__ ":" ASSTR(__LINE__)

/*
 * Macro:[ERROR]
 * print a custom error message and exit the program
 */
#define ERROR(...)                                                                                \
    do                                                                                            \
    {                                                                                             \
        fprintf(stderr, "\n\n*******\n[ERROR](%s)<%s:%d>\n\t", __FUNCTION__, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__);                                                             \
        fprintf(stderr, "\n*******\n");                                                           \
        exit(1);                                                                                  \
    } while (0);

/*
 * Macro:[NEVER]
 * indicates that a point in the code should never be reached
 */
#define NEVER                                                                                 \
    do                                                                                        \
    {                                                                                         \
        ERROR("Control should never reach this point; please report this to the developers.") \
    } while (0);

/*
 * Macro:[ASSERT]
 * evaluate an expression, works the same way as the C-macro assert
 * except that DEBUG does not affect it (it is always active)
 * also prints the file and line info and exits the program
 * if the expression evaluates to false
 */
#define ASSERT(expr)                                                                                                  \
    do                                                                                                                \
    {                                                                                                                 \
        if (!((expr)))                                                                                                \
        {                                                                                                             \
            fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n*******\n", __FILE__, __FUNCTION__, __LINE__, #expr); \
            exit(1);                                                                                                  \
        }                                                                                                             \
    } while (0);

// maps a->0,A->0,c->1,C->1,g->2,G->2,t->3,T=>3,n->4,N->5
int refToInt[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 47
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 79
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 95
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 111
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 127
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 143
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 159
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 175
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 191
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 207
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 223
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 239
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4  // 255
};

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <truth.bcf> <call.bcf>\n", argv[0]);
        return 1;
    }

    htsFile *file1 = bcf_open(argv[1], "r");
    htsFile *file2 = bcf_open(argv[2], "r");

    if (!file1 || !file2)
    {
        fprintf(stderr, "Error opening files.\n");
        return 1;
    }

    bcf_hdr_t *hdr1 = bcf_hdr_read(file1);
    bcf_hdr_t *hdr2 = bcf_hdr_read(file2);
    bcf1_t *rec1 = bcf_init();
    bcf1_t *rec2 = bcf_init();

    int nSamples = bcf_hdr_nsamples(hdr1);
    int *gt_arr1 = NULL, *gt_arr2 = NULL, ngt1 = 0, ngt2 = 0;
    int *nSitesDiscordant = (int *)calloc(nSamples, sizeof(int));

    // number of sites (per individual) where the call file has the site but ind i is missing
    int *nSites_callmis = (int *)calloc(nSamples, sizeof(int));

    int *nSites_compared_forSample = (int *)calloc(nSamples, sizeof(int));
    int f1b1 = 0, f1b2 = 0, f2b1 = 0, f2b2 = 0;
    int f1a1 = 0, f1a2 = 0, f2a1 = 0, f2a2 = 0;
    int sidx1 = 0, sidx2 = 0;

    double missingness_rate = 0.0;
    double discordance_rate = 0.0;
    double concordance_rate = 0.0;
    int nSitesRetained = 0;

    // assume both files have one same contig
    ASSERT(0 == strcmp(bcf_hdr_id2name(hdr1, rec1->rid), bcf_hdr_id2name(hdr2, rec2->rid)));

    int nSites1 = 0;
    int nSites2 = 0;
    int nSitesFile1NotInFile2 = 0;

    while (bcf_read(file1, hdr1, rec1) == 0)
    {
        nSites1++;

        if (bcf_read(file2, hdr2, rec2) != 0)
        {

            // file2 finished
            if (rec2 == NULL)
            {
                // continue with next rec in file1
                continue;
            }
            if (rec1->pos != rec2->pos)
            {
                if (rec1->pos > rec2->pos)
                {
                    // call file cannot contain any sites not included in truth file
                    NEVER;
                }
                else
                {
                    nSitesFile1NotInFile2++;
                    // continue with next rec in file1
                    continue;
                }
                NEVER;
            }
        }
        nSites2++;

        ASSERT(0 == bcf_unpack(rec1, BCF_UN_STR));
        ASSERT(0 == bcf_unpack(rec2, BCF_UN_STR));

        ngt1 = bcf_get_genotypes(hdr1, rec1, &gt_arr1, &ngt1);
        ngt2 = bcf_get_genotypes(hdr2, rec2, &gt_arr2, &ngt2);

        for (int i = 0; i < nSamples; i++)
        {
            sidx1 = 2 * i;
            sidx2 = sidx1 + 1;

            f1a1 = bcf_gt_allele(gt_arr1[sidx1]);
            f1a2 = bcf_gt_allele(gt_arr1[sidx2]);

            f2a1 = bcf_gt_allele(gt_arr2[sidx1]);
            f2a2 = bcf_gt_allele(gt_arr2[sidx2]);

            // true gt file (file1) can never have missing
            ASSERT(f1a1 != -1);
            ASSERT(f1a2 != -1);

            if (-1 == f2a1)
            {
                if (-1 == f2a2)
                {
                    nSites_callmis[i]++;
                    continue;
                }
                else
                {
                    NEVER;
                }
            }

            f1b1 = refToInt[*rec1->d.allele[bcf_gt_allele(gt_arr1[sidx1])]];
            f1b2 = refToInt[*rec1->d.allele[bcf_gt_allele(gt_arr1[sidx2])]];

            f2b1 = refToInt[*rec2->d.allele[bcf_gt_allele(gt_arr2[sidx1])]];
            f2b2 = refToInt[*rec2->d.allele[bcf_gt_allele(gt_arr2[sidx2])]];

            ASSERT(f1b1 >= 0);
            ASSERT(f1b1 <= 3);
            ASSERT(f1b2 >= 0);
            ASSERT(f1b2 <= 3);

            ASSERT(f2b1 >= 0);
            ASSERT(f2b1 <= 3);
            ASSERT(f2b2 >= 0);
            ASSERT(f2b2 <= 3);

            nSites_compared_forSample[i]++;

            if (f1b1 != f2b1)
            {
                if (f1b1 != f2b2)
                {
                    nSitesDiscordant[i]++;
                    continue;
                }
                else
                {
                    // f2b2 used;exclude from check
                    // now check for f1b2
                    if (f1b2 != f2b1)
                    {
                        nSitesDiscordant[i]++;
                        continue;
                    }
                }
            }
            else
            {
                // f2b1 used;exclude from check
                // now check for f1b2
                if (f1b2 != f2b2)
                {
                    nSitesDiscordant[i]++;
                    continue;
                }
            }

        } // samples loop
    }

    fprintf(stderr, "Comparint the genotypes in the truth file %s and the call file %s\n", argv[1], argv[2]);
    fprintf(stderr, "Number of sites in truth file: %d\n", nSites1);
    fprintf(stderr, "Number of sites in the truth file not in the call file: %d\n", nSitesFile1NotInFile2);
    fprintf(stderr, "Number of sites used in comparison: %d\n", nSites2);

    nSitesRetained = nSites1 - nSitesFile1NotInFile2;
    ASSERT(nSitesRetained == nSites2);

    for (int i = 0; i < nSamples; i++)
    {

        // sample, number of sites with discordant genotypes, number of sites compared
        // fprintf(stdout,"%s\t%d\t%d\n", hdr1->samples[i], nSitesDiscordant[i], nSites_compared_forSample[i]);
        fprintf(stdout, "%s", hdr1->samples[i]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d", nSites1);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d", nSitesRetained);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d", nSites_compared_forSample[i]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d", nSites_callmis[i]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d", nSitesDiscordant[i]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d", nSites_compared_forSample[i] - nSitesDiscordant[i]);
        fprintf(stdout, "\t");

        missingness_rate = (double)nSites_compared_forSample[i] / (double)nSites1;
        fprintf(stdout, "%f", missingness_rate);
        fprintf(stdout, "\t");

        discordance_rate = (double)nSitesDiscordant[i] / (double)nSites_compared_forSample[i];
        fprintf(stdout, "%f", discordance_rate);
        fprintf(stdout, "\t");

        concordance_rate = (double)(nSites_compared_forSample[i] - nSitesDiscordant[i]) / (double)nSites_compared_forSample[i];
        fprintf(stdout, "%f", concordance_rate);
        fprintf(stdout, "\n");

        ASSERT(nSites_compared_forSample[i] == nSitesRetained - nSites_callmis[i]);
    }

    free(nSitesDiscordant);
    free(nSites_callmis);
    free(nSites_compared_forSample);
    free(gt_arr1);
    free(gt_arr2);
    bcf_destroy(rec1);
    bcf_destroy(rec2);
    bcf_hdr_destroy(hdr1);
    bcf_hdr_destroy(hdr2);
    bcf_close(file1);
    bcf_close(file2);

    return 0;
}
