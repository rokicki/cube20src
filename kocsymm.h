/*1:*/
#line 21 "./kocsymm.w"

#ifndef KOCSYMM_H
#define KOCSYMM_H
#include "cubepos.h"

/*:1*//*2:*/
#line 79 "./kocsymm.w"

const int CORNERSYMM= 2187;
const int EDGEOSYMM= 2048;
const int EDGEPERM= 495;
/*18:*/
#line 356 "./kocsymm.w"

const int KOCSYMM= 16;
const int CORNERRSYMM= 168;

/*:18*//*19:*/
#line 378 "./kocsymm.w"

struct corner_mapinfo{
unsigned short minbits;
unsigned char csymm,minmap;
};

/*:19*/
#line 83 "./kocsymm.w"

typedef unsigned short lookup_type;
class kocsymm{
public:
kocsymm():csymm(0),eosymm(0),epsymm(0){}
kocsymm(int c,int eo,int ep):csymm(c),eosymm(eo),epsymm(ep){}
/*3:*/
#line 99 "./kocsymm.w"

kocsymm(int):csymm(0),eosymm(0),epsymm(0){init();}
static void init();

/*:3*//*5:*/
#line 115 "./kocsymm.w"

inline bool operator<(const kocsymm&kc)const{
if(csymm!=kc.csymm)return csymm<kc.csymm;
if(eosymm!=kc.eosymm)return eosymm<kc.eosymm;
return epsymm<kc.epsymm;
}
inline bool operator==(const kocsymm&kc)const{
return kc.csymm==csymm&&kc.eosymm==eosymm&&kc.epsymm==epsymm;
}
inline bool operator!=(const kocsymm&kc)const{
return kc.csymm!=csymm||kc.eosymm!=eosymm||kc.epsymm!=epsymm;
}

/*:5*//*9:*/
#line 167 "./kocsymm.w"

void move(int mv){
csymm= cornermove[csymm][mv];
eosymm= edgeomove[eosymm][mv];
epsymm= edgepmove[epsymm][mv];
}

/*:9*//*10:*/
#line 181 "./kocsymm.w"

kocsymm(const cubepos&cp);
void set_coset(cubepos&cp);

/*:10*//*24:*/
#line 468 "./kocsymm.w"

void canon_into(kocsymm&kc)const;

/*:24*//*26:*/
#line 495 "./kocsymm.w"

int calc_symm()const;

/*:26*//*28:*/
#line 518 "./kocsymm.w"

static inline int in_Kociemba_group(int mv){return edgepmove[0][mv]==0;}

/*:28*/
#line 89 "./kocsymm.w"

/*6:*/
#line 134 "./kocsymm.w"

static lookup_type cornermove[CORNERSYMM][NMOVES];
static lookup_type edgeomove[EDGEOSYMM][NMOVES];
static lookup_type edgepmove[EDGEPERM][NMOVES];

/*:6*//*11:*/
#line 208 "./kocsymm.w"

static lookup_type epsymm_compress[1<<12];
static lookup_type epsymm_expand[EDGEOSYMM];

/*:11*//*20:*/
#line 386 "./kocsymm.w"

static lookup_type cornersymm_expand[CORNERRSYMM];
static corner_mapinfo cornersymm[CORNERSYMM];
static lookup_type edgeomap[EDGEOSYMM][KOCSYMM];
static lookup_type edgepmap[EDGEPERM][KOCSYMM];
static lookup_type edgepxor[EDGEPERM][2];

/*:20*/
#line 90 "./kocsymm.w"

lookup_type csymm,eosymm,epsymm;
};

/*:2*//*4:*/
#line 109 "./kocsymm.w"

static kocsymm identity_kc(1);

/*:4*//*29:*/
#line 551 "./kocsymm.w"

const int FACT4= 24;
const int C8_4= 70;
class permcube{
public:
permcube();
/*42:*/
#line 736 "./kocsymm.w"

inline bool operator<(const permcube&pc)const{
if(et!=pc.et)return et<pc.et;
if(em!=pc.em)return em<pc.em;
if(eb!=pc.eb)return eb<pc.eb;
if(etp!=pc.etp)return etp<pc.etp;
if(emp!=pc.emp)return emp<pc.emp;
if(ebp!=pc.ebp)return ebp<pc.ebp;
if(c8_4!=pc.c8_4)return c8_4<pc.c8_4;
if(ctp!=pc.ctp)return ctp<pc.ctp;
return cbp<pc.cbp;
}
inline bool operator==(const permcube&pc)const{
return et==pc.et&&em==pc.em&&eb==pc.eb&&
etp==pc.etp&&emp==pc.emp&&ebp==pc.ebp&&
c8_4==pc.c8_4&&ctp==pc.ctp&&cbp==pc.cbp;
}
inline bool operator!=(const permcube&pc)const{
return et!=pc.et||em!=pc.em||eb!=pc.eb||
etp!=pc.etp||emp!=pc.emp||ebp!=pc.ebp||
c8_4!=pc.c8_4||ctp!=pc.ctp||cbp!=pc.cbp;
}

/*:42*//*45:*/
#line 784 "./kocsymm.w"

void move(int mv);

/*:45*//*47:*/
#line 829 "./kocsymm.w"

void init_edge_from_cp(const cubepos&cp);
void init_corner_from_cp(const cubepos&cp);
permcube(const cubepos&cp);
void set_edge_perm(cubepos&cp)const;
void set_corner_perm(cubepos&cp)const;
void set_perm(cubepos&cp)const;

/*:47*/
#line 557 "./kocsymm.w"

static void init();
/*31:*/
#line 580 "./kocsymm.w"

static unsigned char s4inv[FACT4];
static unsigned char s4mul[FACT4][FACT4];
static unsigned char s4compress[256];
static unsigned char s4expand[FACT4];

/*:31*//*36:*/
#line 668 "./kocsymm.w"

static unsigned char c8_4_compact[256];
static unsigned char c8_4_expand[C8_4];
static unsigned char c8_4_parity[C8_4];

/*:36*//*39:*/
#line 707 "./kocsymm.w"

static unsigned char c12_8[EDGEPERM];
static lookup_type c8_12[C8_4];

/*:39*//*43:*/
#line 772 "./kocsymm.w"

static unsigned short eperm_move[EDGEPERM][NMOVES];
static int cperm_move[C8_4][NMOVES];

/*:43*/
#line 559 "./kocsymm.w"

unsigned short et,em,eb;
unsigned char etp,emp,ebp;
unsigned char c8_4,ctp,cbp;
};

/*:29*//*30:*/
#line 569 "./kocsymm.w"

static permcube identity_pc;

/*:30*//*55:*/
#line 1006 "./kocsymm.w"

#endif

/*:55*/
