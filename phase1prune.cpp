/*2:*/
#line 83 "./phase1prune.w"

#include "phase1prune.h"
#include <iostream> 
#include <cstdio> 
using namespace std;
/*4:*/
#line 112 "./phase1prune.w"

unsigned int phase1prune::memsize;
unsigned char*phase1prune::mem;
int phase1prune::file_checksum;
const char*const phase1prune::filename= "p1p1h.dat";

/*:4*//*18:*/
#line 368 "./phase1prune.w"

static unsigned char map_phase1_offsets[KOCSYMM][3];
static int map_phase1[2][12][256];

/*:18*/
#line 88 "./phase1prune.w"

/*5:*/
#line 126 "./phase1prune.w"

static int datahash(unsigned int*dat,int sz,int seed){
while(sz> 0){
sz-= 4;
seed= 37*seed+*dat++;
}
return seed;
}

/*:5*/
#line 89 "./phase1prune.w"

/*9:*/
#line 163 "./phase1prune.w"

void phase1prune::gen_table(){
memset(mem,-1,memsize);
mem[0]= 0;
int seen= 1;
cout<<"Gen phase1"<<flush;
for(int d= 1;;d++){
int lastiter= (seen==CORNERRSYMM*EDGEOSYMM*EDGEPERM);
int seek= d-1;
int at= 0;
for(int cs= 0;cs<CORNERRSYMM;cs++){
int csymm= kocsymm::cornersymm_expand[cs];
for(int eosymm= 0;eosymm<EDGEOSYMM;eosymm++)
for(int epsymm= 0;epsymm<EDGEPERM;epsymm++,at+= BYTES_PER_ENTRY)
if(mem[at]==seek)
{
/*10:*/
#line 194 "./phase1prune.w"

int deltadist[NMOVES];
for(int mv= 0;mv<NMOVES;mv++){
int rd= 0;
kocsymm kc(csymm,eosymm,epsymm);
kc.move(mv);
corner_mapinfo&cm= kocsymm::cornersymm[kc.csymm];
for(int m= cm.minmap;cm.minbits>>m;m++)
if((cm.minbits>>m)&1){
int deosymm= 
kocsymm::edgeomap[kocsymm::edgepxor[kc.epsymm][m>>3]^kc.eosymm][m];
int depsymm= kocsymm::edgepmap[kc.epsymm][m];
int dat= ((cm.csymm*EDGEOSYMM+deosymm)
*EDGEPERM+depsymm)*BYTES_PER_ENTRY;
rd= mem[dat];
if(rd==255){
rd= d;
mem[dat]= rd;
seen++;
}
}
deltadist[mv]= rd-seek;
}
/*11:*/
#line 232 "./phase1prune.w"

for(int b= 0;b<3;b++){
int v= 0;
int clim= 1;
for(int c= clim;c>=0;c--){
int vv= 0;
int cnts[3];
cnts[0]= cnts[1]= cnts[2]= 0;
for(int t= 2;t>=0;t--){
vv= 2*vv+deltadist[3*b+9*c+t];
cnts[1+deltadist[3*b+9*c+t]]++;
}
if(cnts[0]> 0&&cnts[2]> 0)
error("! bad delta distance values within one face turn set");
if(cnts[0])
vv+= 7;
else
vv+= 8;
v= 16*v+vv;
}
mem[at+b+1]= v;
}

/*:11*/
#line 217 "./phase1prune.w"


/*:10*/
#line 179 "./phase1prune.w"

}
}
cout<<" "<<d<<flush;
if(lastiter)
break;
}
cout<<" done."<<endl<<flush;
}

/*:9*//*12:*/
#line 262 "./phase1prune.w"

const int CHUNKSIZE= 65536;
int phase1prune::read_table(){
FILE*f= fopen(filename,"rb");
if(f==0)
return 0;
int togo= memsize;
unsigned char*p= mem;
int seed= 0;
while(togo> 0){
unsigned int siz= (togo> CHUNKSIZE?CHUNKSIZE:togo);
if(fread(p,1,siz,f)!=siz){
cerr<<"Out of data in "<<filename<<endl;
fclose(f);
return 0;
}
seed= datahash((unsigned int*)p,siz,seed);
togo-= siz;
p+= siz;
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

/*:12*//*13:*/
#line 300 "./phase1prune.w"

void phase1prune::write_table(){
FILE*f= fopen(filename,"wb");
if(f==0)
error("! cannot write pruning file to current directory");
if(fwrite(mem,1,memsize,f)!=memsize)
error("! error writing pruning table");
if(fwrite(&file_checksum,sizeof(int),1,f)!=1)
error("! error writing pruning table");
fclose(f);
}

/*:13*//*14:*/
#line 315 "./phase1prune.w"

void phase1prune::check_integrity(){
if(file_checksum!=datahash((unsigned int*)mem,memsize,0))
error("! integrity of pruning table compromised");
cout<<"Verified integrity of phase one pruning data: "
<<file_checksum<<endl;
}

/*:14*//*17:*/
#line 344 "./phase1prune.w"

int phase1prune::lookup(const kocsymm&kc){
corner_mapinfo&cm= kocsymm::cornersymm[kc.csymm];
int m= cm.minmap;
int r= mem[BYTES_PER_ENTRY*(((cm.csymm*EDGEOSYMM)+
kocsymm::edgeomap[kocsymm::edgepxor[kc.epsymm][m>>3]^kc.eosymm][m])*495+
kocsymm::edgepmap[kc.epsymm][m])];
return r;
}

/*:17*//*19:*/
#line 375 "./phase1prune.w"

int phase1prune::lookup(const kocsymm&kc,int togo,int&nextmovemask){
corner_mapinfo&cm= kocsymm::cornersymm[kc.csymm];
int m= cm.minmap;
int off= 
BYTES_PER_ENTRY*(((cm.csymm*EDGEOSYMM)+
kocsymm::edgeomap[kocsymm::edgepxor[kc.epsymm][m>>3]^kc.eosymm][m])*495+
kocsymm::edgepmap[kc.epsymm][m]);
int r= mem[off];
if(togo<r){
nextmovemask= 0;
}else if(togo> r+1){
nextmovemask= ALLMOVEMASK;
}else{
int(*p)[256]= map_phase1[togo-r];
unsigned char*o= map_phase1_offsets[m];
nextmovemask= p[o[0]][mem[off+1]]+p[o[1]][mem[off+2]]+
p[o[2]][mem[off+3]];
}
return r;
}

/*:19*/
#line 90 "./phase1prune.w"

void phase1prune::init(int suppress_writing){
static int initialized= 0;
if(initialized)
return;
initialized= 1;
/*7:*/
#line 145 "./phase1prune.w"

memsize= BYTES_PER_ENTRY*CORNERRSYMM*EDGEOSYMM*EDGEPERM;
mem= (unsigned char*)malloc(memsize);
if(mem==0)
error("! no memory");

/*:7*//*15:*/
#line 326 "./phase1prune.w"

if(read_table()==0){
gen_table();
file_checksum= datahash((unsigned int*)mem,memsize,0);
if(!suppress_writing)
write_table();
}

/*:15*//*20:*/
#line 403 "./phase1prune.w"

for(int m= 0;m<KOCSYMM;m++){
for(int f= 0;f<3;f++){
int mv= f*TWISTS;
int mv2= cubepos::move_map[m][mv];
int f2= mv2/TWISTS;
int key= 0;
if(mv2%TWISTS==TWISTS-1)
key++;
if(f2>=3)
key+= 2;
key+= 4*(f2%3);
map_phase1_offsets[cubepos::invm[m]][f]= key;
}
}

/*:20*//*21:*/
#line 423 "./phase1prune.w"

for(int slack= 0;slack<2;slack++){
for(int key= 0;key<12;key++){
int nv[16];
for(int nyb= 0;nyb<16;nyb++){
int bits= 0;
if(slack&&nyb<=7){
bits= 7;
}else if(slack==0&&nyb>=7){
bits= 0;
}else{
bits= 7-(nyb&7);
}
if(key&1)
bits= ((bits&1)<<2)+(bits&2)+(bits>>2);
if(key&2)
bits<<= 3*TWISTS;
bits<<= TWISTS*(key>>2);
nv[nyb]= bits;
}
int*a= map_phase1[slack][key];
for(int byte= 0;byte<256;byte++)
a[byte]= nv[byte&15]|(((nv[byte>>4]<<(3*TWISTS))|
(nv[byte>>4]>>(3*TWISTS)))&0777777);
}
}

/*:21*/
#line 96 "./phase1prune.w"

}

/*:2*//*22:*/
#line 452 "./phase1prune.w"

moveseq phase1prune::solve(kocsymm kc){
moveseq r;
int d= phase1prune::lookup(kc);
while(d> 0){
int nmm= 0;
int t= phase1prune::lookup(kc,d,nmm);
if(t==0)
break;
if(t!=d)
error("! did not make progress");
if(nmm==0)
error("! no solution?");
int mv= ffs(nmm)-1;
r.push_back(mv);
kc.move(mv);
d--;
}
return r;
}

/*:22*/
