/*3:*/
#line 167 "./phase2prune.w"

#include "phase2prune.h"
#include <iostream> 
#include <cstdio> 
using namespace std;
/*4:*/
#line 202 "./phase2prune.w"

struct corner_reduce{
unsigned char m,parity;
lookup_type c,minbits;
}corner_reduction[FACT8];
lookup_type edgeud_remap[KOCSYMM][FACT8];

/*:4*//*7:*/
#line 233 "./phase2prune.w"

int phase2prune::cornermax;
unsigned int phase2prune::memsize;
unsigned int*phase2prune::mem;

/*:7*//*18:*/
#line 459 "./phase2prune.w"

const char*const phase2prune::filename= "p2p1h.dat";
int phase2prune::file_checksum;

/*:18*/
#line 172 "./phase2prune.w"

/*5:*/
#line 214 "./phase2prune.w"

inline int corner_coordinate(const permcube&pc){
return(pc.c8_4*FACT4+pc.ctp)*FACT4+pc.cbp;
}
inline int edge_coordinate(const permcube&pc){
return(permcube::c12_8[pc.et]*FACT4+pc.etp)*FACT4+pc.ebp;
}

/*:5*//*19:*/
#line 467 "./phase2prune.w"

static int datahash(unsigned int*dat,int sz,int seed){
while(sz> 0){
sz-= 4;
seed= 37*seed+*dat++;
}
return seed;
}

/*:19*/
#line 173 "./phase2prune.w"

/*11:*/
#line 314 "./phase2prune.w"

int phase2prune::lookup(const cubepos&cp){
permcube pc(cp);
return lookup(pc);
}
int phase2prune::lookup(const permcube&pc){
int cc= corner_coordinate(pc);
corner_reduce&cr= corner_reduction[cc];
int off= cr.c*FACT8+edgeud_remap[cr.m][edge_coordinate(pc)];
int r= (mem[off>>3]>>(4*(off&7)))&0xf;
if(r==0&&pc==identity_pc)
return 0;
else
return r+1;
}

/*:11*//*13:*/
#line 353 "./phase2prune.w"

void phase2prune::gen_table(){
memset(mem,255,memsize);
cout<<"Gen phase2"<<flush;
mem[0]&= ~14;
int seen= 1;
for(int d= 0;d<15;d++){
unsigned int seek= (d?d-1:1);
int newval= d;
for(int c8_4= 0;c8_4<C8_4;c8_4++)
for(int ctp= 0;ctp<FACT4;ctp++)
for(int cbp= 0;cbp<FACT4;cbp++){
permcube pc;
pc.c8_4= c8_4;
pc.ctp= ctp;
pc.cbp= cbp;
int oc= corner_coordinate(pc);
corner_reduce&cr= corner_reduction[oc];
if(cr.minbits&1){
/*14:*/
#line 390 "./phase2prune.w"

permcube pc2,pc3,pc4;
cubepos cp2,cp3;
int off= corner_reduction[oc].c*(FACT8/8);
for(int mv= 0;mv<NMOVES;mv++){
if(!kocsymm::in_Kociemba_group(mv))
continue;
pc2= pc;
pc2.move(mv);
int dest_off= corner_coordinate(pc2);
corner_reduce&cr= corner_reduction[dest_off];
int destat= cr.c*(FACT8/8);
for(int m= cr.m;(1<<m)<=cr.minbits;m++)
if((cr.minbits>>m)&1){
/*15:*/
#line 413 "./phase2prune.w"

int at= 0;
for(int e8_4= 0;e8_4<C8_4;e8_4++){
int et= permcube::c8_12[e8_4];
int t1= permcube::eperm_move[et][mv];
int eb= kocsymm::epsymm_compress[0xf0f-kocsymm::epsymm_expand[et]];
int t2= permcube::eperm_move[eb][mv]&31;
int dst1= permcube::c12_8[t1>>5]*24*24;
t1&= 31;
for(int etp= 0;etp<FACT4;etp++)
for(int ebp= 0;ebp<FACT4;ebp++,at++){
if(mem[off+(at>>3)]==0xffffffff){
ebp+= 7;
at+= 7;
}else if(((mem[off+(at>>3)]>>(4*(at&7)))&0xf)==seek){
/*16:*/
#line 436 "./phase2prune.w"

int etp1= permcube::s4mul[etp][t1];
int ebp1= permcube::s4mul[ebp][t2];
int dat= edgeud_remap[m][dst1+etp1*24+ebp1];
int val= (mem[destat+(dat>>3)]>>(4*(dat&7)))&0xf;
if(val==0xf){
mem[destat+(dat>>3)]-= (0xf-newval)<<(4*(dat&7));
seen++;
}

/*:16*/
#line 428 "./phase2prune.w"

}
}
}

/*:15*/
#line 404 "./phase2prune.w"

}
}

/*:14*/
#line 372 "./phase2prune.w"
;
}
}
#ifndef QUARTER
if(d==0)
mem[0]&= ~15;
#endif
cout<<" "<<d<<flush;
}
cout<<" done."<<endl<<flush;
}

/*:13*//*20:*/
#line 482 "./phase2prune.w"

const int CHUNKSIZE= 65536;
int phase2prune::read_table(){
FILE*f= fopen(filename,"rb");
if(f==0)
return 0;
int togo= memsize;
unsigned int*p= mem;
int seed= 0;
while(togo> 0){
unsigned int siz= (togo> CHUNKSIZE?CHUNKSIZE:togo);
if(fread(p,1,siz,f)!=siz){
cerr<<"Out of data in "<<filename<<endl;
fclose(f);
return 0;
}
seed= datahash(p,siz,seed);
togo-= siz;
p+= siz/sizeof(unsigned int);
}
if(fread(&file_checksum,sizeof(int),1,f)!=1){
cerr<<"Out of data in "<<filename<<endl;
fclose(f);
return 0;
}
fclose(f);
if(file_checksum!=seed){
cerr<<"Bad checksum in "<<filename<<"; expected "
<<file_checksum<<" but saw "<<seed<<endl;
return 0;
}
return 1;
}

/*:20*//*21:*/
#line 520 "./phase2prune.w"

void phase2prune::write_table(){
FILE*f= fopen(filename,"wb");
if(f==0)
error("! cannot write pruning file to current directory");
if(fwrite(mem,1,memsize,f)!=memsize)
error("! error writing pruning table");
if(fwrite(&file_checksum,sizeof(int),1,f)!=1)
error("! error writing pruning table");
fclose(f);
}

/*:21*//*22:*/
#line 535 "./phase2prune.w"

void phase2prune::check_integrity(){
if(file_checksum!=datahash(mem,memsize,0))
error("! integrity of pruning table compromised");
cout<<"Verified integrity of phase two pruning data: "
<<file_checksum<<endl;
}

/*:22*//*25:*/
#line 572 "./phase2prune.w"

moveseq phase2prune::solve(const permcube&pc,int maxlen){
moveseq r;
for(int d= lookup(pc);d<=maxlen;d++)
if(solve(pc,d,CANONSEQSTART,r)){
reverse(r.begin(),r.end());
break;
}
return r;
}
int phase2prune::solve(const permcube&pc,int togo,
int canonstate,moveseq&r){
if(lookup(pc)> togo)
return 0;
if(pc==identity_pc)
return 1;
if(togo--<=0)
return 0;
permcube pc2;
int mask= cubepos::cs_mask(canonstate)&0227227227;
while(mask){
int ntogo= togo;
int mv= ffs(mask)-1;
mask&= mask-1;
pc2= pc;
pc2.move(mv);
if(solve(pc2,ntogo,cubepos::next_cs(canonstate,mv),r)){
r.push_back(mv);
return 1;
}
}
return 0;
}


/*:25*/
#line 174 "./phase2prune.w"

void phase2prune::init(int suppress_writing){
static int initialized= 0;
if(initialized)
return;
initialized= 1;
/*8:*/
#line 240 "./phase2prune.w"

cornermax= 0;
for(int c8_4= 0;c8_4<C8_4;c8_4++)
for(int ctp= 0;ctp<FACT4;ctp++)
for(int cbp= 0;cbp<FACT4;cbp++){
permcube pc;
pc.c8_4= c8_4;
pc.ctp= ctp;
pc.cbp= cbp;
int oc= corner_coordinate(pc);
int minc= oc;
int minm= 0;
int minbits= 1;
cubepos cp;
pc.set_perm(cp);
for(int m= 1;m<16;m++){
cubepos cp2;
cp.remap_into(m,cp2);
permcube pc2(cp2);
int tc= corner_coordinate(pc2);
if(tc<minc){
minc= tc;
minm= m;
minbits= 1<<m;
}else if(tc==minc)
minbits|= 1<<m;
}
corner_reduce&cr= corner_reduction[oc];
if(oc==minc)
cr.c= cornermax++;
cr.m= minm;
cr.c= corner_reduction[minc].c;
cr.minbits= minbits;
cr.parity= (permcube::c8_4_parity[c8_4]+ctp+cbp)&1;
}

/*:8*//*9:*/
#line 277 "./phase2prune.w"

int at= 0;
cubepos cp,cp2;
for(int e8_4= 0;e8_4<C8_4;e8_4++){
permcube pc;
pc.et= permcube::c8_12[e8_4];
pc.eb= kocsymm::epsymm_compress[0xf0f-kocsymm::epsymm_expand[pc.et]];
for(int etp= 0;etp<FACT4;etp++){
pc.etp= etp;
for(int ebp= 0;ebp<FACT4;ebp++,at++){
pc.ebp= ebp;
for(int m= 0;m<KOCSYMM;m++){
pc.set_edge_perm(cp);
cp.remap_into(m,cp2);
permcube pc2(cp2);
edgeud_remap[m][at]= edge_coordinate(pc2);
}
}
}
}

/*:9*//*10:*/
#line 301 "./phase2prune.w"

memsize= cornermax*(FACT8/2);
mem= (unsigned int*)malloc(memsize);
if(mem==0)
error("! no memory in phase2prune");

/*:10*//*23:*/
#line 546 "./phase2prune.w"

if(read_table()==0){
gen_table();
file_checksum= datahash(mem,memsize,0);
if(!suppress_writing)
write_table();
}

/*:23*/
#line 180 "./phase2prune.w"

}

/*:3*/
