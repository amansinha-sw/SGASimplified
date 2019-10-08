#ifndef READ_H
#define READ_H
typedef struct{
    char *seq;
    char *qual;
    int len;
}read;

const char alphabet[] = "ACGT";
#endif
