#include <stdio.h>
#include <stdlib.h>
#include "read.h"
#include "fm_index.h"
#include "util.h"

#define m_minSupportLowQuality 3
#define m_minSupportHighQuality 2
#define m_highQualityCutoff 20

char* get_rev_comp(char *seq, int length){
    char *rev = (char*)malloc(length*sizeof(int));
    for(int i=0;i<length;i++){
        rev[length-1-i] = seq[i];
        if(rev[length-1-i]=='A')
            rev[length-1-i] = 'T';
        else if(rev[length-1-i]=='C')
            rev[length-1-i] = 'G';
        else if(rev[length-1-i]=='G')
            rev[length-1-i] = 'C';
        else if(rev[length-1-i]=='T')
            rev[length-1-i] = 'A';
    }
    printf("%s<->%s\n", seq, rev);
    return rev;
}

int kmer_correct(int i, int k_idx, read *item, int threshold, fm_index *index, int kmer_length){
    char untrusted_base = item->seq[i];
    char trusted_base = '$';
    char kmer[kmer_length];
    int base_idx = i - k_idx;
    for(int j=0;j<kmer_length;j++)
        kmer[j] = item->seq[j+k_idx];
    int max_count = 0;
    for(int j=0;j<4;j++){
        if(alphabet[j]!=untrusted_base)
            kmer[base_idx] = alphabet[j];
        else
            continue;
        int count = get_kmer_count(kmer, kmer_length, index);
        count += get_kmer_count(get_rev_comp(kmer, kmer_length), kmer_length, index);
        if(count>=threshold){
            if(trusted_base == '$'){
                trusted_base = alphabet[j];
                max_count = count;
            }
            else{
                printf("MUltiple values\n");
                return 0;
            }
        }
    }
    if(trusted_base != '$' && max_count >= threshold){
        item->seq[i] = trusted_base;
        printf("FOund correct base for position %d:%c threshold=%d count=%d : %s using kmer = %s\n",i,trusted_base, threshold,max_count, item->seq, kmer);
        return 1;
    }
    return 0;
}

void read_correct(read *item, fm_index *index, int kmer_length, int max_attempts){
    int *trust_value = (int*)calloc(item->len,sizeof(int));
    int num_kmers = item->len - kmer_length + 1;
    int *min_phred = (int*)malloc(num_kmers*sizeof(int));
    int all_trusted = 1;
    int done =0;
    int round = 0;
    while(done!=1 && round<max_attempts){
        round++;
        //char *rev_comp = get_rev_comp(item->seq, item->len);
        for(int i=0;i<num_kmers;i++){
            int min_phred = 512;
            for(int j=0;j<kmer_length;j++){
                min_phred = min_phred > item->qual[j+i] ? item->qual[j+i] : min_phred;
            }
            int threshold = m_minSupportLowQuality;
            if(min_phred >= m_highQualityCutoff)
                threshold = m_minSupportHighQuality;
            int count = get_kmer_count(&item->seq[i], kmer_length, index);
            printf("kmer count for %s : %d\n",&item->seq[i], count);
            count += get_kmer_count(get_rev_comp(&item->seq[i], kmer_length), kmer_length, index);
            printf("rev kmer count for %s : %d\n",get_rev_comp(&item->seq[i], kmer_length), count);
            printf("kmer=%d count=%d threshold=%d\n", i, count, threshold);
            if(count >= threshold){
                for(int j=0;j<kmer_length;j++){
                    trust_value[i+j] = 1;
                }
            }
        }
        all_trusted = 1;
        for(int i=0;i<item->len;i++)
            all_trusted *= trust_value[i];

        int corrected = 0;
        if(all_trusted!=1){
            for(int i=0;i<item->len;i++){
                printf("%c: Trusted=%d\n",item->seq[i],trust_value[i]);
            
                if(trust_value[i]!=1){
                    int phred = item->qual[i];
                    int threshold = m_minSupportLowQuality;
                    if(phred >= m_highQualityCutoff)
                        threshold = m_minSupportHighQuality;
                    int left_kmer_k = i + 1 >= kmer_length ? i + 1 - kmer_length : 0;
                    corrected = kmer_correct(i, left_kmer_k, item, threshold, index, kmer_length);
                    if(corrected==1)
                        break;
                    int right_kmer_k = i < item->len - kmer_length ? i : item->len - kmer_length;
                    corrected = kmer_correct(i, right_kmer_k, item, threshold, index, kmer_length);
                    if(corrected==1)
                        break;
                }    
            }
            if(corrected!=1){
                printf("No more correction possible on this read\n");
                done = 1;
                break;
            }
            else
                continue;
        }
        else
            break;
    }
    if(all_trusted == 1)
        printf("Corrected read : %s\n",item->seq);
}

int main(int argc, char *argv[]){
    fm_index index;
    //printf("Reading fm-index\n");
    init_fm_index_from_file(&index, argv[1]);
    read *reads = init_reads_from_file(argv[2]);
    read_correct(&reads[0], &index, 3, 5);
}
