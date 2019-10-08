#include <stdio.h>
#include <stdlib.h>
#include "fm_index.h"

//#define PRINT_FM_INDEX 0

void init_fm_index_from_file(fm_index *index, char *filename){
    FILE *fp;
    fp = fopen(filename, "r");
    int length;
    fscanf(fp, "%d\n", &length);
    fscanf(fp, "%d %d %d %d %d %d\n",&index->alphabet_starts[0],&index->alphabet_starts[1],&index->alphabet_starts[2],&index->alphabet_starts[3],&index->alphabet_starts[4],&index->alphabet_starts[5]);
    index->len = length;
    index->l_column = (char*)malloc(length*sizeof(char));
    fread(index->l_column, 1, length, fp);
    index->ranks.len = length;
    index->ranks.stride = 1;
    for(int i=0;i<5;i++)
        index->ranks.counts[i]=(int*)calloc(length,sizeof(int));
    int j = -1;
    char c = index->l_column[0];
    if(c=='A')
        j = 1;
    else if(c=='C')
        j = 2;
    else if(c=='G')
        j = 3;
    else if(c=='T')
        j = 4;
    else if(c=='$')
        j = 0;
    index->ranks.counts[j][0] = 1;
    #ifdef PRINT_FM_INDEX
    printf("0(%c): %d %d %d %d %d\n",c,index->ranks.counts[0][0], index->ranks.counts[1][0], index->ranks.counts[2][0], index->ranks.counts[3][0], index->ranks.counts[4][0]);
    #endif
    for(int i=1;i<length;i++){
        for(int k=0;k<5;k++)
            index->ranks.counts[k][i] = index->ranks.counts[k][i-1];
        c = index->l_column[i];
        if(c=='A')
            j = 1;
        else if(c=='C')
            j = 2;
        else if(c=='G')
            j = 3;
        else if(c=='T')
            j = 4;
        else if(c=='$')
            j = 0;
        index->ranks.counts[j][i] = index->ranks.counts[j][i] + 1;
        #ifdef PRINT_FM_INDEX
        printf("%d(%c): %d %d %d %d %d\n",i,c,index->ranks.counts[0][i], index->ranks.counts[1][i], index->ranks.counts[2][i], index->ranks.counts[3][i], index->ranks.counts[4][i]);
        #endif
    }
}

read* init_reads_from_file(char *filename){
    int reads_count;
    FILE *fp;
    fp = fopen(filename, "r");
    fscanf(fp, "%d\n", &reads_count);
    read *reads = (read*) malloc(reads_count*sizeof(read));
    for(int i=0;i<reads_count;i++){
        fscanf(fp, "%d\n", &reads[i].len);
        reads[i].seq = (char*)malloc((reads[i].len+1)*sizeof(char));
        reads[i].qual = (char*)malloc((reads[i].len+1)*sizeof(char));
        fread(reads[i].seq, 1, reads[i].len+1, fp);
        fread(reads[i].qual, 1, reads[i].len+1, fp);
        reads[i].seq[reads[i].len] = '\0';
        reads[i].qual[reads[i].len] = '\0';
        printf("len =%d : %s %s\n",reads[i].len,reads[i].seq, reads[i].qual);
    }
    return reads;
}
