/*
    Full-featured high-performance C implementation for:
    1. SNP overlap filtering
    2. CpG-disrupting SNP filtering
    3. InDel flank filtering
    4. Depth filtering
    5. Context-aware methylation filtering
    Compile:
        gcc -O3 -march=native taps_filter.c -o taps_filter
    Run:
        ./taps_filter methylation.tsv variants.vcf output.tsv
    Optional:
        ./taps_filter_full methylation.tsv variants.vcf output.tsv 5 5

        flank=5
        min_depth=5
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_LINE 8192
#define HASH_SIZE 2000003

typedef struct SNPNode {

    char chrom[64];
    int pos;

    char ref;
    char alt;

    struct SNPNode* next;

} SNPNode;

typedef struct {

    int start;
    int end;

} Interval;

typedef struct {

    Interval* data;

    int size;
    int capacity;

} IntervalList;


/* =========================
   GLOBALS
========================= */
SNPNode* snp_hash[HASH_SIZE];

typedef struct ChromIntervals {

    char chrom[64];

    IntervalList intervals;

    struct ChromIntervals* next;

} ChromIntervals;

ChromIntervals* interval_map = NULL;

/* =========================
   HASH FUNCTIONS
========================= */
unsigned hash_key(const char* chrom, int pos) {

    unsigned h = (unsigned)pos;

    while (*chrom) {

        h = h * 131 + *chrom;
        chrom++;
    }

    return h % HASH_SIZE;
}

/* =========================
   SNP FUNCTIONS
========================= */
void add_snp(
    const char* chrom,
    int pos,
    char ref,
    char alt
) {

    unsigned h = hash_key(chrom, pos);

    SNPNode* node = malloc(sizeof(SNPNode));

    strcpy(node->chrom, chrom);

    node->pos = pos;

    node->ref = ref;
    node->alt = alt;

    node->next = snp_hash[h];

    snp_hash[h] = node;
}

SNPNode* get_snp(
    const char* chrom,
    int pos
) {

    unsigned h = hash_key(chrom, pos);

    SNPNode* node = snp_hash[h];

    while (node) {

        if (
            node->pos == pos &&
            strcmp(node->chrom, chrom) == 0
        ) {
            return node;
        }

        node = node->next;
    }

    return NULL;
}


/* =========================
   INTERVAL FUNCTIONS
========================= */

ChromIntervals* get_chrom_intervals(
    const char* chrom
) {

    ChromIntervals* cur = interval_map;

    while (cur) {

        if (strcmp(cur->chrom, chrom) == 0)
            return cur;

        cur = cur->next;
    }

    ChromIntervals* node =
        malloc(sizeof(ChromIntervals));

    strcpy(node->chrom, chrom);

    node->intervals.size = 0;
    node->intervals.capacity = 1024;

    node->intervals.data =
        malloc(sizeof(Interval) * node->intervals.capacity);

    node->next = interval_map;

    interval_map = node;

    return node;
}

void add_interval(
    const char* chrom,
    int start,
    int end
) {

    ChromIntervals* ci =
        get_chrom_intervals(chrom);

    IntervalList* list = &ci->intervals;

    if (list->size >= list->capacity) {

        list->capacity *= 2;

        list->data = realloc(
            list->data,
            sizeof(Interval) * list->capacity
        );
    }

    list->data[list->size].start = start;
    list->data[list->size].end = end;

    list->size++;
}

int cmp_interval(
    const void* a,
    const void* b
) {

    Interval* x = (Interval*)a;
    Interval* y = (Interval*)b;

    return x->start - y->start;
}

void sort_intervals() {

    ChromIntervals* cur = interval_map;

    while (cur) {

        qsort(
            cur->intervals.data,
            cur->intervals.size,
            sizeof(Interval),
            cmp_interval
        );

        cur = cur->next;
    }
}


/* =========================
   BINARY SEARCH OVERLAP
========================= */

int overlap_indel(
    const char* chrom,
    int pos
) {

    ChromIntervals* ci = interval_map;

    while (ci) {

        if (strcmp(ci->chrom, chrom) == 0)
            break;

        ci = ci->next;
    }

    if (!ci)
        return 0;

    IntervalList* list = &ci->intervals;

    int left = 0;
    int right = list->size - 1;

    while (left <= right) {

        int mid = (left + right) / 2;

        Interval* iv = &list->data[mid];

        if (pos < iv->start) {

            right = mid - 1;

        } else if (pos > iv->end) {

            left = mid + 1;

        } else {

            return 1;
        }
    }

    return 0;
}


/* =========================
   LOAD VCF
========================= */

void load_vcf(
    const char* vcf_file,
    int flank
) {

    FILE* fp = fopen(vcf_file, "r");

    if (!fp) {

        fprintf(stderr,
            "ERROR: cannot open VCF\n");

        exit(1);
    }

    char line[MAX_LINE];

    while (fgets(line, MAX_LINE, fp)) {

        if (line[0] == '#')
            continue;

        char* fields[20];

        int i = 0;

        char* token = strtok(line, "\t");

        while (token) {

            fields[i++] = token;
            token = strtok(NULL, "\t");
        }

        char* chrom = fields[0];

        int pos = atoi(fields[1]);

        char* ref = fields[3];
        char* alt = fields[4];

        double qual = atof(fields[5]);

        if (qual < 30)
            continue;

        /* DP parse */

        char* format = fields[8];
        char* sample = fields[9];

        int dp = 0;

        char fmt_copy[1024];
        char sample_copy[1024];

        strcpy(fmt_copy, format);
        strcpy(sample_copy, sample);

        char* fmt_tok = strtok(fmt_copy, ":");
        char* sample_tok = strtok(sample_copy, ":");

        while (fmt_tok && sample_tok) {

            if (strcmp(fmt_tok, "DP") == 0) {

                dp = atoi(sample_tok);
                break;
            }

            fmt_tok = strtok(NULL, ":");
            sample_tok = strtok(NULL, ":");
        }

        if (dp < 5)
            continue;

        /* SNP */

        if (
            strlen(ref) == 1 &&
            strlen(alt) == 1
        ) {

            add_snp(
                chrom,
                pos,
                ref[0],
                alt[0]
            );
        }

        /* InDel */

        else {

            int start = pos - flank;
            int end =
                pos +
                (int)(
                    strlen(ref) > strlen(alt)
                    ? strlen(ref)
                    : strlen(alt)
                ) +
                flank;

            add_interval(
                chrom,
                start,
                end
            );
        }
    }

    fclose(fp);

    sort_intervals();
}


/* =========================
   CONTEXT FILTER
========================= */

int should_remove_by_snp(
    const char* chrom,
    int pos,
    const char* context
) {

    /* direct overlap */

    if (get_snp(chrom, pos))
        return 1;

    /* CpG-disrupting */

    if (
        context[0] == 'C' &&
        context[1] == 'G'
    ) {

        SNPNode* partner =
            get_snp(chrom, pos + 1);

        if (
            partner &&
            partner->ref == 'G' &&
            partner->alt == 'A'
        ) {
            return 1;
        }
    }

    return 0;
}


/* =========================
   FILTER METHYLATION
========================= */

void filter_methylation(
    const char* meth_file,
    const char* out_file,
    int min_depth
) {

    FILE* fin = fopen(meth_file, "r");

    FILE* fout = fopen(out_file, "w");

    if (!fin || !fout) {

        fprintf(stderr,
            "ERROR: cannot open files\n");

        exit(1);
    }

    char line[MAX_LINE];

    long total = 0;
    long retained = 0;

    /* header */

    fgets(line, MAX_LINE, fin);

    fprintf(fout, "%s", line);

    while (fgets(line, MAX_LINE, fin)) {

        total++;

        char buffer[MAX_LINE];

        strcpy(buffer, line);

        char* fields[20];

        int i = 0;

        char* token = strtok(buffer, "\t");

        while (token) {

            fields[i++] = token;
            token = strtok(NULL, "\t");
        }

        char* chrom = fields[0];

        int pos = atoi(fields[1]) + 1;

        char* context = fields[8];

        int depth = atoi(fields[11]);

        if (depth < min_depth)
            continue;

        if (
            should_remove_by_snp(
                chrom,
                pos,
                context
            )
        ) {
            continue;
        }

        if (
            overlap_indel(
                chrom,
                pos
            )
        ) {
            continue;
        }

        fprintf(fout, "%s", line);

        retained++;
    }

    fclose(fin);
    fclose(fout);

    fprintf(stderr,
        "\nTotal: %ld\nRetained: %ld\n",
        total,
        retained
    );
}


/* =========================
   MAIN
========================= */

int main(
    int argc,
    char* argv[]
) {

    if (argc < 4) {

        fprintf(stderr,
            "\nUsage:\n"
            "%s methylation.tsv variants.vcf output.tsv [flank] [min_depth]\n\n",
            argv[0]
        );

        return 1;
    }

    int flank = 5;
    int min_depth = 5;

    if (argc >= 5)
        flank = atoi(argv[4]);

    if (argc >= 6)
        min_depth = atoi(argv[5]);

    fprintf(stderr,
        "\nLoading variants...\n");

    load_vcf(
        argv[2],
        flank
    );

    fprintf(stderr,
        "Filtering methylation...\n");

    filter_methylation(
        argv[1],
        argv[3],
        min_depth
    );

    fprintf(stderr,
        "Done.\n");

    return 0;
}
