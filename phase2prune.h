/*1:*/
#line 138 "./phase2prune.w"

#ifndef PHASE2PRUNE_H
#define PHASE2PRUNE_H
#include "kocsymm.h"

/*:1*//*2:*/
#line 150 "./phase2prune.w"

const int FACT8= 40320;
class phase2prune{
public:
static void init(int suppress_writing= 0);
static int lookup(const cubepos&cp);
static int lookup(const permcube&pc);
/*12:*/
#line 340 "./phase2prune.w"

static void gen_table();
static int read_table();
static void write_table();
static void check_integrity();

/*:12*//*24:*/
#line 560 "./phase2prune.w"

static moveseq solve(const permcube&pc,int maxlen= 30);
static moveseq solve(const cubepos&cp,int maxlen= 30){
permcube pc(cp);
return solve(pc,maxlen);
}
static int solve(const permcube&pc,int togo,int canonstate,moveseq&seq);

/*:24*/
#line 157 "./phase2prune.w"

/*6:*/
#line 226 "./phase2prune.w"

static int cornermax;
static unsigned int memsize;
static unsigned int*mem;

/*:6*//*17:*/
#line 452 "./phase2prune.w"

static const char*const filename;
static int file_checksum;

/*:17*/
#line 158 "./phase2prune.w"

};
#endif

/*:2*/
