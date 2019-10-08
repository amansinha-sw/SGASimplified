#ifndef FM_INDEX_H
#define FM_INDEX_H

typedef struct{
    int *counts[5];//$,A,C,G,T
    int len;
    int stride;
}alphabet_ranks;

typedef struct{
    char *l_column;
    int len;
    alphabet_ranks ranks;
    int alphabet_starts[6];//$,A,C,G,T,-
}fm_index;

void init_interval(int *interval_start, int *interval_end, char c, fm_index *index){
    int i = -1;
    if(c=='A')
        i = 1;
    else if(c=='C')
        i = 2;
    else if(c=='G')
        i = 3;
    else if(c=='T')
        i = 4;
    *interval_start = index->alphabet_starts[i];
    *interval_end = index->alphabet_starts[i+1]-1;
}

void update_interval(int *interval_start, int *interval_end, char c, fm_index *index){
    int i = -1;
    if(c=='A')
        i = 1;
    else if(c=='C')
        i = 2;
    else if(c=='G')
        i = 3;
    else if(c=='T')
        i = 4;
    
    *interval_start = index->alphabet_starts[i] + index->ranks.counts[i][*interval_start-1];
    *interval_end = index->alphabet_starts[i] + index->ranks.counts[i][*interval_end]-1;
}

int get_kmer_count(char *kmer, int length, fm_index *index){
    int kmer_start = 0, kmer_end = 0;
    init_interval(&kmer_start, &kmer_end, kmer[length-1], index);
    for(int i=length-2; i>=0; i--)
        update_interval(&kmer_start, &kmer_end, kmer[i], index);
    return kmer_end - kmer_start + 1;
}
#endif
