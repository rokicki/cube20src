/*1:*/
#line 66 "./phase1prune.w"

#ifndef PHASE1PRUNE_H
#define PHASE1PRUNE_H
#include "kocsymm.h"
/*6:*/
#line 137 "./phase1prune.w"

const int BYTES_PER_ENTRY= 4;

/*:6*/
#line 70 "./phase1prune.w"
;
class phase1prune{
public:
static void init(int suppress_writing= 0);
static int lookup(const kocsymm&kc,int&mask);
/*8:*/
#line 154 "./phase1prune.w"

static void gen_table();
static int read_table();
static void write_table();
static void check_integrity();

/*:8*//*16:*/
#line 336 "./phase1prune.w"

static int lookup(const kocsymm&kc);
static int lookup(const kocsymm&kc,int togo,int&nextmovemask);
static moveseq solve(kocsymm kc);

/*:16*/
#line 75 "./phase1prune.w"

/*3:*/
#line 103 "./phase1prune.w"

static unsigned int memsize;
static unsigned char*mem;
static int file_checksum;
static const char*const filename;

/*:3*/
#line 76 "./phase1prune.w"

};
#endif

/*:1*/
