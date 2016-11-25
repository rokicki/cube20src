/*1:*/
#line 17 "./hcoset.w"

const char*BANNER= 
"This is hcoset 1.0, (C) 2010 Tomas Rokicki.  All Rights Reserved.";
#include "phase1prune.h"
#include <pthread.h> 
#include <iostream> 
#include <map> 
#include <set> 
using namespace std;
/*2:*/
#line 48 "./hcoset.w"

int verbose= 1;
int numthreads= 1;
const int MAX_THREADS= 32;

/*:2*//*4:*/
#line 85 "./hcoset.w"

int skipwrite= 0;

/*:4*//*6:*/
#line 103 "./hcoset.w"

moveseq repseq;
set<permcube> world;
kocsymm repkc;
permcube reppc;
cubepos repcp;

/*:6*//*10:*/
#line 157 "./hcoset.w"

int slow= 2;
int maxsearchdepth= 35;
int maxdepth= 35;
int global_depth;

/*:10*//*15:*/
#line 281 "./hcoset.w"

const int FACT8= 40320;
const int PAGESIZE= (FACT8*FACT4/2/8);
unsigned char**bitp1,**bitp2;

/*:15*//*17:*/
#line 299 "./hcoset.w"

long long uniq= 0;
long long probes= 0;
const long long TARGET= FACT8*(long long)FACT8*(long long)FACT4/2;
#ifdef LEVELCOUNTS
long long uniq_ulev= 0;
long long sum_ulev[30];
#endif

/*:17*//*19:*/
#line 343 "./hcoset.w"

#ifdef FASTCLEAN
unsigned char touched[FACT8];
int did_a_prepass;
#endif

/*:19*//*21:*/
#line 370 "./hcoset.w"

unsigned char permtobit[FACT4];
unsigned char bittoperm[FACT4];

/*:21*//*23:*/
#line 390 "./hcoset.w"

unsigned char saveb;
unsigned char*savep;

/*:23*//*29:*/
#line 548 "./hcoset.w"

int disable_prepass= 0;

/*:29*//*32:*/
#line 594 "./hcoset.w"

const int SQMOVES= 3;
short rearrange[2][SQMOVES][1<<12];

/*:32*//*35:*/
#line 643 "./hcoset.w"

const int PREPASS_MOVES= 10;
unsigned short eperm_map[FACT8/2][PREPASS_MOVES];

/*:35*//*39:*/
#line 780 "./hcoset.w"

#ifdef LEVELCOUNTS
unsigned char bc[1<<12];
#endif

/*:39*//*43:*/
#line 841 "./hcoset.w"

const int STRIDE= 16;

/*:43*//*44:*/
#line 852 "./hcoset.w"

#include "corner_order.h"

/*:44*//*49:*/
#line 965 "./hcoset.w"

pthread_mutex_t mutex;

/*:49*//*56:*/
#line 1076 "./hcoset.w"

unsigned char use_count[FACT8];
int work_done;

/*:56*//*61:*/
#line 1190 "./hcoset.w"

int this_level_did_prepass= 0;

/*:61*//*66:*/
#line 1292 "./hcoset.w"

int search_terminated_early= 0;
int dont_count_after_max_search= 0;
int need_count_bits= 0;
long long enoughbits= TARGET;
int fast20;

/*:66*//*75:*/
#line 1600 "./hcoset.w"

#include "bestsol.h"

/*:75*//*76:*/
#line 1613 "./hcoset.w"

int first_coset;
const int TOTALCOSETS= 55882296;
int coset_count;

/*:76*//*79:*/
#line 1665 "./hcoset.w"

map<int,int> symcount;

/*:79*//*82:*/
#line 1804 "./hcoset.w"

#ifdef LEVELCOUNTS
long long levprev,levuniq,levsum;
char levmul[FACT8];
int k3map[FACT8];
#endif

/*:82*/
#line 26 "./hcoset.w"

/*13:*/
#line 208 "./hcoset.w"

void slowsearch1(const kocsymm&kc,const permcube&pc,int togo,
int movemask,int canon){
if(togo==0){
if(kc==identity_kc){
probes++;
world.insert(pc);
}
return;
}
togo--;
kocsymm kc2;
permcube pc2;
int newmovemask;
while(movemask){
int mv= ffs(movemask)-1;
movemask&= movemask-1;
kc2= kc;
kc2.move(mv);
int nd= phase1prune::lookup(kc2,togo,newmovemask);
if(nd<=togo){
pc2= pc;
pc2.move(mv);
int new_canon= cubepos::next_cs(canon,mv);
slowsearch1(kc2,pc2,togo,
newmovemask&cubepos::cs_mask(new_canon),
new_canon);
}
}
}

/*:13*//*14:*/
#line 243 "./hcoset.w"

void slowsearch1(const kocsymm&kc,const permcube&pc){
duration();
for(int d= phase1prune::lookup(repkc);d<maxsearchdepth;d++){
probes= 0;
long long prevlev= uniq;
slowsearch1(kc,pc,d,ALLMOVEMASK,CANONSEQSTART);
uniq= world.size();
long long thislev= uniq-prevlev;
if(verbose)
cout<<"Tests at "<<d<<" "<<probes<<" in "<<duration()
<<" uniq "<<uniq<<" lev "<<thislev<<endl;
}
}

/*:14*//*18:*/
#line 314 "./hcoset.w"

vector<unsigned char*> pageq;
unsigned char*getpage(){
unsigned char*r= 0;
if(pageq.size()> 0){
r= pageq[pageq.size()-1];
pageq.pop_back();
}else{
r= (unsigned char*)malloc(PAGESIZE+1);
if(r==0)
error("! no memory");
}
return r;
}
unsigned char*getclearedpage(){
unsigned char*r= getpage();
memset(r,0,PAGESIZE);
return r;
}
void freepage(unsigned char*r){
pageq.push_back(r);
}

/*:18*//*24:*/
#line 402 "./hcoset.w"

void flushbit(){
if(savep!=0){
if(0==(*savep&saveb)){
*savep|= saveb;
uniq++;
}
savep= 0;
}
}
void setonebit(const permcube&pc){
int cindex= (pc.c8_4*FACT4+pc.ctp)*FACT4+pc.cbp;
unsigned int eindex= (((permcube::c12_8[pc.et]*FACT4)+pc.etp)*FACT4/2+
(pc.ebp>>1))*FACT4+permtobit[pc.emp];
#ifdef FASTCLEAN
touched[cindex]= 1;
#endif
probes++;
flushbit();
savep= bitp1[cindex]+(eindex>>3);
__builtin_prefetch(savep);
saveb= 1<<(eindex&7);
}

/*:24*//*25:*/
#line 434 "./hcoset.w"

void slowsearch2(const kocsymm&kc,const permcube&pc,int togo,
int movemask,int canon){
if(togo==0){
if(kc==identity_kc)
setonebit(pc);
return;
}
togo--;
kocsymm kc2;
permcube pc2;
int newmovemask;
while(movemask){
int mv= ffs(movemask)-1;
movemask&= movemask-1;
kc2= kc;
kc2.move(mv);
int nd= phase1prune::lookup(kc2,togo,newmovemask);
if(nd<=togo){
pc2= pc;
pc2.move(mv);
int new_canon= cubepos::next_cs(canon,mv);
int movemask3= newmovemask&cubepos::cs_mask(new_canon);
if(togo==1){
permcube pc3;
while(movemask3){
int mv2= ffs(movemask3)-1;
movemask3&= movemask3-1;
pc3= pc2;
pc3.move(mv2);
setonebit(pc3);
}
}else{
slowsearch2(kc2,pc2,togo,movemask3,new_canon);
}
}
}
}

/*:25*//*27:*/
#line 481 "./hcoset.w"

void slowsearch2(const kocsymm&kc,const permcube&pc){
duration();
for(int d= phase1prune::lookup(repkc);d<maxsearchdepth;d++){
probes= 0;
long long prevlev= uniq;
slowsearch2(kc,pc,d,ALLMOVEMASK,CANONSEQSTART);
flushbit();
long long thislev= uniq-prevlev;
if(verbose)
cout<<"Tests at "<<d<<" "<<probes<<" in "<<duration()
<<" uniq "<<uniq<<" lev "<<thislev<<endl;
}
}

/*:27*//*28:*/
#line 499 "./hcoset.w"

void unpack_edgecoord(permcube&pc,int e8_4,int epp1,int epp2){
pc.et= permcube::c8_12[e8_4];
pc.etp= epp1;
pc.ebp= epp2;
pc.eb= kocsymm::epsymm_compress[0xf0f-kocsymm::epsymm_expand[pc.et]];
pc.em= 0;
}
void unpack_edgecoord(permcube&pc,int coord){
unpack_edgecoord(pc,coord/(FACT4*FACT4),coord/FACT4%FACT4,coord%FACT4);
}

/*:28*//*37:*/
#line 712 "./hcoset.w"

void innerloop3(unsigned char*dst,unsigned char**srcs,int base){
dst+= 3*base;
unsigned char*end= dst+12*24*3;
unsigned char tval= *end;
unsigned short*cpp= eperm_map[base];
for(;dst<end;dst+= 3,cpp+= PREPASS_MOVES){
int wf2= *(int*)(srcs[3]+cpp[3]);
int wr2= *(int*)(srcs[4]+cpp[4]);
int wb2= *(int*)(srcs[8]+cpp[8]);
int wl2= *(int*)(srcs[9]+cpp[9]);
*(int*)dst= 
(((wf2&0xfff)|
rearrange[0][0][wr2&0xfff]|
rearrange[0][1][wb2&0xfff]|
rearrange[0][2][wl2&0xfff])<<12)|
((wf2>>12)&0xfff)|
rearrange[1][0][(wr2>>12)&0xfff]|
rearrange[0][1][(wb2>>12)&0xfff]|
rearrange[1][2][(wl2>>12)&0xfff]|
*(int*)(srcs[0]+cpp[0])|
*(int*)(srcs[1]+cpp[1])|
*(int*)(srcs[2]+cpp[2])|
*(int*)(srcs[5]+cpp[5])|
*(int*)(srcs[6]+cpp[6])|
*(int*)(srcs[7]+cpp[7]);
}
*end= tval;
}

/*:37*//*38:*/
#line 746 "./hcoset.w"

int countbits(int*a){
int r= 0;
const int mask1= 0x55555555;
const int mask2= 0x33333333;
const int mask3= 0x0f0f0f0f;
for(int i= 0;i<PAGESIZE;i+= 24){
int w1= *a++;
int w2= *a++;
int w3= *a++;
w1= (w1&mask1)+((w1>>1)&mask1)+(w2&mask1);
w2= ((w2>>1)&mask1)+(w3&mask1)+((w3>>1)&mask1);
int s1= (w1&mask2)+((w1>>2)&mask2)+
(w2&mask2)+((w2>>2)&mask2);
s1= (s1&mask3)+((s1>>4)&mask3);
w1= *a++;
w2= *a++;
w3= *a++;
w1= (w1&mask1)+((w1>>1)&mask1)+(w2&mask1);
w2= ((w2>>1)&mask1)+(w3&mask1)+((w3>>1)&mask1);
int s2= (w1&mask2)+((w1>>2)&mask2)+
(w2&mask2)+((w2>>2)&mask2);
s1+= (s2&mask3)+((s2>>4)&mask3);
r+= 255&((s1>>24)+(s1>>16)+(s1>>8)+s1);
}
return r;
}

/*:38*//*41:*/
#line 797 "./hcoset.w"

#ifdef LEVELCOUNTS
int parity(int coord){
return(permcube::c8_4_parity[coord/(FACT4*FACT4)]^(coord/24)^coord)&1;
}
int countbits2(int cperm,int*a,int&rv2){
int coparity= parity(cperm);
int r= 0,r2= 0;
int ind= 0;
for(int e8_4= 0;e8_4<C8_4;e8_4++){
int p2= coparity^permcube::c8_4_parity[e8_4];
for(int epp1= 0;epp1<FACT4;epp1++){
int off1= (p2^epp1)&1;
int off2= 1-off1;
for(int epp2= 0;epp2<FACT4;epp2+= 2,ind+= 2){
int w= *a;
int v1= bc[(w&0xfff)];
int v2= bc[((w>>12)&0xfff)];
r+= v1+v2;
r2+= v1*levmul[ind+off1]+v2*levmul[ind+off2];
a= (int*)(((char*)a)+3);
}
}
}
rv2= r2;
return r;
}
#endif

/*:41*//*42:*/
#line 832 "./hcoset.w"

struct elemdata{
unsigned char*dst;
unsigned char*from[PREPASS_MOVES];
unsigned char e84map[70];
};

/*:42*//*45:*/
#line 857 "./hcoset.w"

void calcneighbors(int cperm,int*a){
permcube pc,pc2;
pc.c8_4= cperm/(FACT4*FACT4);
pc.ctp= cperm/FACT4%FACT4;
pc.cbp= cperm%FACT4;
for(int mv= 0;mv<NMOVES;mv++){
if(!kocsymm::in_Kociemba_group(mv))
continue;
pc2= pc;
pc2.move(mv);
*a++= (pc2.c8_4*FACT4+pc2.ctp)*FACT4+pc2.cbp;
}
}

/*:45*//*46:*/
#line 877 "./hcoset.w"

int moveseq16[]= {2,9,0,9,2,9,0,0,0,11,2,11,0,11,2,0};
unsigned char initorder[70];
void doouter(int r){
elemdata edata[16];
int neighbors[PREPASS_MOVES];
int tcperm= cornerorder[r];
permcube pc,pc2;
for(int i= 0;i<STRIDE;i++){
elemdata*e= edata+i;
int cperm= cornerorder[r+i];
if(cperm!=tcperm)
error("! inconsistent corner order");
calcneighbors(cperm,neighbors);
int dp= 0;
for(int mv= 0;mv<NMOVES;mv++){
if(!kocsymm::in_Kociemba_group(mv))
continue;
e->from[dp]= bitp2[neighbors[dp]];
if(mv==moveseq16[i])
tcperm= neighbors[dp];
dp++;
}
if(i==0){
for(int j= 0;j<70;j++)
e->e84map[j]= initorder[j];
}else{
pc= identity_pc;
for(int j= 0;j<70;j++){
pc.c8_4= (e-1)->e84map[j];
pc.move(moveseq16[i-1]);
e->e84map[j]= pc.c8_4;
}
}
e->dst= bitp1[cperm];
}
for(int i= 0;i<70;i++)
for(int j= 0;j<STRIDE;j++){
elemdata*e= edata+j;
innerloop3(e->dst,e->from,e->e84map[i]*12*24);
}
}

/*:46*//*51:*/
#line 975 "./hcoset.w"

void get_global_lock(){
pthread_mutex_lock(&mutex);
}
void release_global_lock(){
pthread_mutex_unlock(&mutex);
}

/*:51*//*52:*/
#line 988 "./hcoset.w"

FILE*singfile;
int singcount;
int singlevel= 99;
void showsing(const permcube&pc){
cubepos cp,cp2;
pc.set_perm(cp);
cubepos::mul(cp,repcp,cp2);
if(singfile)
fprintf(singfile,"SING %s\n",cp2.Singmaster_string());
else
cout<<"SING "<<cp2.Singmaster_string()<<endl;
singcount++;
}

/*:52*//*55:*/
#line 1036 "./hcoset.w"

void showunset(int cperm){
int*pbits= (int*)bitp1[cperm];
permcube pc;
pc.c8_4= cperm/(FACT4*FACT4);
pc.ctp= cperm/FACT4%FACT4;
pc.cbp= cperm%FACT4;
int coparity= (permcube::c8_4_parity[pc.c8_4]^pc.ctp^pc.cbp)&1;
for(int i= 0;i<FACT8/8*3;i++){
if(pbits[i]!=-1){
int t= ~pbits[i];
while(t){
int j= ffs(t)-1;
t&= t-1;
int k= (i<<5)+j;
int ep8_4= k/(24*24*12);
int epp1= k/(24*12)%24;
int epp2= k/24%12*2;
int emperm= bittoperm[k%24];
int par0= (permcube::c8_4_parity[ep8_4]^epp1^
emperm^coparity)&1;
unpack_edgecoord(pc,ep8_4,epp1,epp2+par0);
pc.emp= emperm;
get_global_lock();
showsing(pc);
release_global_lock();
}
}
}
}

/*:55*//*57:*/
#line 1087 "./hcoset.w"

int get_prepass_work(){
get_global_lock();
int r= work_done;
if(r<FACT8){
for(int i= 0;i<STRIDE;i++)
bitp1[cornerorder[r+i]]= getpage();
work_done+= STRIDE;
}
release_global_lock();
if(r>=FACT8)
return-1;
return r;
}

/*:57*//*58:*/
#line 1106 "./hcoset.w"

void finish_prepass_work(int cperm){
int thisblock= 0;
#ifdef LEVELCOUNTS
int this_ulev= 0;
#endif
if(need_count_bits){
#ifdef LEVELCOUNTS
if(need_count_bits> 1)
thisblock= countbits2(cperm,(int*)bitp1[cperm],this_ulev);
else
#endif
thisblock= countbits((int*)bitp1[cperm]);
}
if(global_depth+1>=singlevel&&maxsearchdepth<global_depth)
showunset(cperm);
int neighbors[PREPASS_MOVES];
calcneighbors(cperm,neighbors);
get_global_lock();
uniq+= thisblock;
#ifdef LEVELCOUNTS
uniq_ulev+= this_ulev;
#endif
for(int i= 0;i<PREPASS_MOVES;i++){
if((--use_count[neighbors[i]])==0){
freepage(bitp2[neighbors[i]]);
bitp2[neighbors[i]]= 0;
}
}
release_global_lock();
}

/*:58*//*59:*/
#line 1141 "./hcoset.w"

void*do_prepass_work(void*){
while(1){
int r= get_prepass_work();
if(r<0)
break;
doouter(r);
for(int i= 0;i<STRIDE;i++)
finish_prepass_work(cornerorder[r+i]);
}
return 0;
}

/*:59*//*60:*/
#line 1159 "./hcoset.w"

void doprepass(){
uniq= 0;
#ifdef FASTCLEAN
did_a_prepass= 1;
#endif
#ifdef LEVELCOUNTS
uniq_ulev= 0;
#endif
swap(bitp1,bitp2);
memset(use_count,PREPASS_MOVES,sizeof(use_count));
work_done= 0;
pthread_t p_thread[MAX_THREADS];
for(int ti= 1;ti<numthreads;ti++)
pthread_create(&(p_thread[ti]),NULL,do_prepass_work,0);
do_prepass_work(0);
for(int ti= 1;ti<numthreads;ti++)
pthread_join(p_thread[ti],0);
if(need_count_bits==0)
cout<<"Prepass at "<<global_depth<<" done in "<<duration()
<<"; unique ?"<<endl<<flush;
else
cout<<"Prepass at "<<global_depth<<" done in "<<duration()
<<"; unique "<<uniq<<endl<<flush;
}

/*:60*//*62:*/
#line 1200 "./hcoset.w"

const int CHECKABITSIZE= 64;
struct checkabit{
unsigned char*p;
#ifdef LEVELCOUNTS
unsigned char weight;
#endif
unsigned char b;
};

/*:62*//*69:*/
#line 1380 "./hcoset.w"

int search_work_seq;
int get_search_work(){
get_global_lock();
int r= -1;
while(1){
r= search_work_seq++;
if(r>=NMOVES*NMOVES){
r= -1;
break;
}
int mv1= r/NMOVES;
int mv2= r%NMOVES;
int s= cubepos::next_cs(CANONSEQSTART,mv1);
int mask= cubepos::cs_mask(s);
if((mask>>mv2)&1)
break;
}
release_global_lock();
return r;
}

/*:69*//*78:*/
#line 1643 "./hcoset.w"

const int U2= 1;
void docoset(int seq,const char*movestring);
void docoverelement(int seq,const kocsymm&kc){
int d= phase1prune::lookup(kc);
if(d> maxdepth)
return;
moveseq moves= phase1prune::solve(kc);
moves= cubepos::invert_sequence(moves);
if(moves.size()==0){
moves.push_back(U2);
moves.push_back(U2);
}
char buf[160];
strcpy(buf,cubepos::moveseq_string(moves));
docoset(seq,buf);
}

/*:78*//*80:*/
#line 1671 "./hcoset.w"

int bestsollookup[EDGEOSYMM*EDGEPERM];
unsigned int orderkc(const kocsymm&kc){
return(kc.epsymm<<11)+kc.eosymm;
}
int lookupkc(const kocsymm&kc){
return bestsollookup[(kc.epsymm<<11)+kc.eosymm];
}
void genseqs(int lo,int hi){
phase1prune::init();
cubepos cp,cp2;
int tot= 0;
int c= 0;
for(int eo= 0;eo<EDGEOSYMM;eo++)
for(int ep= 0;ep<EDGEPERM;ep++){
kocsymm kc(0,eo,ep);
kocsymm kc2;
kc.canon_into(kc2);
if(kc==kc2)
bestsollookup[orderkc(kc)]= bestsol[c++];
else
bestsollookup[orderkc(kc)]= bestsollookup[orderkc(kc2)];
}
for(int match= 2;match> 0;match--){
int c= 0;
for(int eo= 0;eo<EDGEOSYMM;eo++){
for(int ep= 0;ep<EDGEPERM;ep++){
kocsymm kc(0,eo,ep);
kocsymm kc2;
kc.canon_into(kc2);
if(!(kc==kc2))
continue;
if(bestsol[c]!=match){
c++;
continue;
}
int cnt= 1;
kc.set_coset(cp);
for(int m= 1;m<16;m++){
cp.remap_into(m,cp2);
kocsymm kc2(cp2);
if(kc2==kc)
cnt|= 1<<m;
}
if(tot+CORNERSYMM<lo||tot>=hi){
if(cnt==1){
tot+= CORNERSYMM;
c++;
continue;
}else if(symcount.find(cnt)!=symcount.end()){
tot+= symcount[cnt];
c++;
continue;
}
}
int tcnt= 0;
for(int co= 0;co<CORNERSYMM;co++){
kc.csymm= co;
int okay= 1;
kc.set_coset(cp);
for(int m= 1;okay&&m<16;m++){
if(0==((cnt>>m)&1))
continue;
cp.remap_into(m,cp2);
kocsymm kc2(cp2);
if(kc2<kc)
okay= 0;
}
if(okay){
if(tot>=lo&&tot<hi)
docoverelement(tot,kc);
tcnt++;
tot++;
}
}
symcount[cnt]= tcnt;
c++;
}
}
}
if(tot!=TOTALCOSETS)
error("! mistake in computation of total cosets");
}

/*:80*//*83:*/
#line 1813 "./hcoset.w"

#ifdef LEVELCOUNTS
void setupk3map(){
int lookaside[1<<12];
memset(lookaside,-1,sizeof(lookaside));
int ind= 0;
cubepos cp;
permcube pc;
for(int e8_4= 0;e8_4<C8_4;e8_4++){
for(int epp1= 0;epp1<FACT4;epp1++){
for(int epp2= 0;epp2<FACT4;epp2++,ind++){
unpack_edgecoord(pc,e8_4,epp1,epp2);
pc.set_perm(cp);
int key= (1<<(cp.e[0]/2))|(1<<(cp.e[3]/2))|
(1<<(cp.e[8]/2))|(1<<(cp.e[11]/2));
if(lookaside[key]<0)
lookaside[key]= ind;
k3map[ind]= lookaside[key];
}
}
}
}
#endif

/*:83*//*85:*/
#line 1847 "./hcoset.w"

#ifdef LEVELCOUNTS
void setup_levmul(const kocsymm&kc,const moveseq&moves){
int x= kc.calc_symm();
memset(levmul,48/x,sizeof(levmul));
int ind= 0;
for(int e8_4= 0;e8_4<70;e8_4++)
for(int epp1= 0;epp1<FACT4;epp1++)
for(int epp2= 0;epp2<FACT4;epp2++,ind++){
permcube pc;
if(k3map[ind]!=ind){
levmul[ind]= levmul[k3map[ind]];
}else{
cubepos cpt,cp2,cp3;
unpack_edgecoord(pc,e8_4,epp1,epp2);
pc.set_perm(cp2);
for(unsigned int i= 0;i<moves.size();i++)
cp2.move(moves[i]);
for(int i= 0;i<8;i++)
cp2.c[i]= cubepos::corner_val(i,0);
cp2.invert_into(cpt);
cpt.remap_into(16,cp2);
kocsymm kc3(cp2);
kocsymm kc1;
kocsymm kc2(cpt);
kc2.canon_into(kc1);
kc3.canon_into(kc2);
cpt.remap_into(32,cp2);
kocsymm kct(cp2);
kct.canon_into(kc3);
if((orderkc(kc2)<orderkc(kc1)&&lookupkc(kc2))||
(orderkc(kc3)<orderkc(kc1)&&lookupkc(kc3))){
levmul[ind]= 0;
}else{
int d= 1;
if(kc1==kc2)
d++;
if(kc1==kc3)
d++;
levmul[ind]= 48/(x*d);
}
}
}
}
#endif

/*:85*/
#line 27 "./hcoset.w"

/*48:*/
#line 956 "./hcoset.w"

struct worker_thread{
/*63:*/
#line 1212 "./hcoset.w"

long long local_probes;
#ifdef LEVELCOUNTS
long long local_ulev;
#endif
int bitcount;
checkabit q[CHECKABITSIZE];
void local_bitflush(){
get_global_lock();
int this_uniq= 0;
#ifdef LEVELCOUNTS
int this_ulev= 0;
#endif
for(int i= 0;i<bitcount;i++){
checkabit&c= q[i];
if(0==(*c.p&c.b)){
*c.p|= c.b;
#ifdef LEVELCOUNTS
this_ulev+= c.weight;
#endif
this_uniq++;
}
}
uniq+= this_uniq;
if(uniq> enoughbits&&global_depth==maxsearchdepth)
search_terminated_early= 1;
#ifdef LEVELCOUNTS
local_ulev+= this_ulev;
#endif
release_global_lock();
local_probes+= bitcount;
bitcount= 0;
}
void local_initialize(){
local_probes= 0;
#ifdef LEVELCOUNTS
local_ulev= 0;
#endif
bitcount= 0;
/*65:*/
#line 1279 "./hcoset.w"

#ifdef FASTCLEAN
memset(local_touched,0,sizeof(local_touched));
#endif

/*:65*/
#line 1251 "./hcoset.w"
;
}

/*:63*//*64:*/
#line 1256 "./hcoset.w"

unsigned char local_touched[FACT8];
void local_setonebit(const permcube&pc){
int cindex= (pc.c8_4*FACT4+pc.ctp)*FACT4+pc.cbp;
int epindex= ((permcube::c12_8[pc.et]*FACT4)+pc.etp)*FACT4+pc.ebp;
int eindex= (epindex>>1)*FACT4+permtobit[pc.emp];
#ifdef FASTCLEAN
local_touched[cindex]= 1;
#endif
if(bitcount>=CHECKABITSIZE)
local_bitflush();
checkabit&c= q[bitcount];
c.p= bitp1[cindex]+(eindex>>3);
__builtin_prefetch(c.p);
#ifdef LEVELCOUNTS
c.weight= levmul[epindex];
#endif
c.b= 1<<(eindex&7);
bitcount++;
}

/*:64*//*68:*/
#line 1333 "./hcoset.w"

void search(const kocsymm&kc,const permcube&pc,int togo,
int movemask,int canon){
if(togo==0){
if(kc==identity_kc)
local_setonebit(pc);
return;
}
togo--;
kocsymm kc2;
permcube pc2;
int newmovemask;
if(search_terminated_early)
return;
while(movemask){
int mv= ffs(movemask)-1;
movemask&= movemask-1;
kc2= kc;
kc2.move(mv);
int nd= phase1prune::lookup(kc2,togo,newmovemask);
if(nd<=togo&&(
togo==nd||togo+nd>=5||
!this_level_did_prepass)){
pc2= pc;
pc2.move(mv);
int new_canon= cubepos::next_cs(canon,mv);
int movemask3= newmovemask&cubepos::cs_mask(new_canon);
if(togo==1){
permcube pc3;
while(movemask3){
int mv2= ffs(movemask3)-1;
movemask3&= movemask3-1;
pc3= pc2;
pc3.move(mv2);
local_setonebit(pc3);
}
}else{
search(kc2,pc2,togo,movemask3,new_canon);
}
}
}
}

/*:68*//*70:*/
#line 1404 "./hcoset.w"

void dowork(){
local_initialize();
if(global_depth<=3){
search(repkc,reppc,global_depth,ALLMOVEMASK,CANONSEQSTART);
}else{
while(1){
int movepair= get_search_work();
if(movepair<0)
break;
int mv1= movepair/NMOVES;
int mv2= movepair%NMOVES;
kocsymm kc(repkc);
permcube pc(reppc);
kc.move(mv1);
pc.move(mv1);
kc.move(mv2);
pc.move(mv2);
int s= cubepos::next_cs(CANONSEQSTART,mv1);
s= cubepos::next_cs(s,mv2);
int mask;
phase1prune::lookup(kc,global_depth-2,mask);
search(kc,pc,global_depth-2,mask&cubepos::cs_mask(s),s);
}
}
local_bitflush();
get_global_lock();
probes+= local_probes;
#ifdef LEVELCOUNTS
uniq_ulev+= local_ulev;
#endif
#ifdef FASTCLEAN
for(int i= 0;i<FACT8;i++)
touched[i]|= local_touched[i];
#endif
release_global_lock();
}

/*:70*/
#line 958 "./hcoset.w"
;
char pad[128];
}workers[MAX_THREADS];

/*:48*//*72:*/
#line 1451 "./hcoset.w"

void*do_search_work(void*s){
worker_thread*w= (worker_thread*)s;
w->dowork();
return 0;
}

/*:72*//*73:*/
#line 1464 "./hcoset.w"

void search(){
duration();
for(int d= phase1prune::lookup(repkc);d<=maxdepth;d++){
global_depth= d;
probes= 0;
long long prevlev= uniq;
#ifdef LEVELCOUNTS
long long prev_ulev= uniq_ulev;
#endif
this_level_did_prepass= !disable_prepass&&d> 1&&
(uniq> 6000000||d> maxsearchdepth);
if(this_level_did_prepass){
need_count_bits= 
(d<=maxsearchdepth||!dont_count_after_max_search);
#ifdef LEVELCOUNTS
if(d<=maxsearchdepth)
need_count_bits++;
#endif
doprepass();
}else{
need_count_bits= 1;
}
if(uniq==TARGET)
break;
search_work_seq= 0;
search_terminated_early= 0;
int did_full_search= 0;
if(fast20&&d==16)
enoughbits= 167000000+uniq/3;
if(d<=maxsearchdepth){
did_full_search= 1;
if(d<=3){
workers[0].dowork();
}else{
pthread_t p_thread[MAX_THREADS];
for(int ti= 1;ti<numthreads;ti++)
pthread_create(&(p_thread[ti]),NULL,do_search_work,
&workers[ti]);
workers[0].dowork();
for(int ti= 1;ti<numthreads;ti++)
pthread_join(p_thread[ti],0);
}
}
if(search_terminated_early){
cout<<"Terminated search at level "<<d
<<" early due to enoughbits"<<endl;
if(dont_count_after_max_search)
need_count_bits= 0;
did_full_search= 0;
}
long long thislev= uniq-prevlev;
#ifdef LEVELCOUNTS
long long thisulev= uniq_ulev-prev_ulev;
#endif
if(global_depth+1>=singlevel&&maxsearchdepth>=global_depth)
for(int cperm= 0;cperm<FACT8;cperm++)
showunset(cperm);
if(verbose){
#ifdef LEVELCOUNTS
if(need_count_bits==0)
cout<<"Tests at "<<d<<" "<<probes<<" in "<<duration()
<<" uniq ? lev ?"<<endl;
else if(did_full_search)
cout<<"Tests at "<<d<<" "<<probes<<" in "<<duration()
<<" uniq "<<uniq<<" lev "<<thislev<<" utot "
<<uniq_ulev<<" thisulev "<<thisulev<<endl;
else
cout<<"Tests at "<<d<<" "<<probes<<" in "<<duration()
<<" uniq "<<uniq<<" lev "<<thislev<<endl;
#else
if(did_full_search)
cout<<"Tests at "<<d<<" "<<probes<<" in "<<duration()
<<" uniq ?"<<endl;
else
cout<<"Tests at "<<d<<" "<<probes<<" in "<<duration()
<<" uniq "<<uniq<<" lev "<<thislev<<endl;
#endif
}
#ifdef LEVELCOUNTS
if(did_full_search)
sum_ulev[d]+= thisulev;
#endif
}
}

/*:73*/
#line 28 "./hcoset.w"

void docoset(int seq,const char*movestring){
/*8:*/
#line 123 "./hcoset.w"

int oldsingcount= singcount;
const char*tmp= movestring;
repseq= cubepos::parse_moveseq(tmp);
if(*tmp)
error("! extra stuff at end of input moveseq");
cout<<"Coset representative "<<seq<<" "<<movestring<<endl;
phase1prune::init(skipwrite);
double cosetstart= walltime();

/*:8*//*9:*/
#line 135 "./hcoset.w"

repkc= identity_kc;
reppc= identity_pc;
repcp= identity_cube;
for(unsigned int i= 0;i<repseq.size();i++){
repkc.move(repseq[i]);
reppc.move(repseq[i]);
repcp.move(repseq[i]);
}
#ifdef LEVELCOUNTS
setup_levmul(repkc,repseq);
#endif

/*:9*//*12:*/
#line 201 "./hcoset.w"

if(slow==0)
slowsearch1(repkc,reppc);

/*:12*//*20:*/
#line 353 "./hcoset.w"

for(int i= 0;i<FACT8;i++)
if(bitp1[i]==0)
bitp1[i]= getclearedpage();
uniq= 0;
#ifdef FASTCLEAN
memset(touched,0,sizeof(touched));
did_a_prepass= 0;
#endif
#ifdef LEVELCOUNTS
uniq_ulev= 0;
#endif

/*:20*//*26:*/
#line 475 "./hcoset.w"

if(slow==1)
slowsearch2(repkc,reppc);

/*:26*//*71:*/
#line 1445 "./hcoset.w"

if(slow> 1)
search();

/*:71*//*74:*/
#line 1554 "./hcoset.w"

int delta= singcount-oldsingcount;
if(singfile)
cout<<"Wrote "<<delta<<" sing positions."<<endl;
for(int i= 0;i<FACT8;i++){
#ifdef FASTCLEAN
if(bitp1[i]!=0&&(touched[i]||did_a_prepass)){
freepage(bitp1[i]);
bitp1[i]= 0;
}
#else
if(bitp1[i]!=0){
freepage(bitp1[i]);
bitp1[i]= 0;
}
#endif
if(bitp2[i]!=0){
freepage(bitp2[i]);
bitp2[i]= 0;
}
}
cout<<"Finished in "<<(walltime()-cosetstart)<<endl;

/*:74*/
#line 30 "./hcoset.w"

}
int main(int argc,char*argv[]){
double progstart= walltime();
duration();
/*3:*/
#line 55 "./hcoset.w"

int oargc= argc;
char**oargv= argv;
while(argc> 1&&argv[1][0]=='-'){
argc--;
argv++;
switch(argv[0][1]){
case'v':verbose++;break;
case'q':verbose= 0;break;
case't':
if(argc<2)
error("! not enough arguments to -t");
if(sscanf(argv[1],"%d",&numthreads)!=1)
error("! bad thread count argument");
if(numthreads<1||numthreads> MAX_THREADS)
error("! bad value for thread count");
argc--;
argv++;
break;
/*5:*/
#line 90 "./hcoset.w"

case'W':skipwrite++;break;

/*:5*//*11:*/
#line 165 "./hcoset.w"

case's':
if(argc<2)
error("! not enough arguments to -s");
if(sscanf(argv[1],"%d",&slow)!=1)
error("! bad -s argument");
if(slow<0||slow> 2)
error("! bad value for -s");
argc--;
argv++;
break;
case'd':
if(argc<2)
error("! not enough arguments to -d");
if(sscanf(argv[1],"%d",&maxdepth)!=1)
error("! bad -d argument");
if(maxdepth<0)
error("! bad value for -d");
if(maxdepth<maxsearchdepth)
maxsearchdepth= maxdepth;
argc--;
argv++;
break;
case'S':
if(argc<2)
error("! not enough arguments to -S");
if(sscanf(argv[1],"%d",&maxsearchdepth)!=1)
error("! bad -S argument");
if(maxsearchdepth<0)
error("! bad value for -S");
argc--;
argv++;
break;

/*:11*//*30:*/
#line 553 "./hcoset.w"

case'U':disable_prepass++;break;

/*:30*//*53:*/
#line 1005 "./hcoset.w"

case'f':
if(argc<2)
error("! not enough arguments to -f");
singfile= fopen(argv[1],"w");
if(singfile==0)
error("! can't open singfile");
argc--;
argv++;
break;
case'L':
if(argc<2)
error("! not enough arguments to -L");
if(sscanf(argv[1],"%d",&singlevel)!=1)
error("! non-numeric argument to -L");
argc--;
argv++;
break;

/*:53*//*67:*/
#line 1301 "./hcoset.w"

case'e':
if(argc<2)
error("! not enough arguments to -e");
if(sscanf(argv[1],"%lld",&enoughbits)!=1)
error("! bad arguments to -e");
argc--;
argv++;
break;
case'a':
dont_count_after_max_search++;
break;
case'F':
fast20++;
dont_count_after_max_search++;
maxsearchdepth= 16;
maxdepth= 20;
break;

/*:67*//*77:*/
#line 1620 "./hcoset.w"

case'r':
if(argc<3)
error("! not enough arguments to -r");
if(sscanf(argv[1],"%d",&first_coset)!=1||
sscanf(argv[2],"%d",&coset_count)!=1)
error("! bad arguments to -r");
argc-= 2;
argv+= 2;
if(first_coset<0||first_coset>=TOTALCOSETS)
error("! bad first value to -r");
if(coset_count<=0)
error("! bad coset count specified to -r");
if(coset_count+first_coset> TOTALCOSETS)
coset_count= TOTALCOSETS-first_coset;
break;

/*:77*/
#line 74 "./hcoset.w"

default:
error("! bad argument");
}
}

/*:3*/
#line 35 "./hcoset.w"

/*7:*/
#line 112 "./hcoset.w"

if(verbose)
cout<<BANNER<<endl<<flush;
for(int i= 0;i<oargc;i++)
cout<<" "<<oargv[i];
cout<<endl;

/*:7*//*16:*/
#line 288 "./hcoset.w"

bitp1= (unsigned char**)calloc(FACT8,sizeof(unsigned char*));
bitp2= (unsigned char**)calloc(FACT8,sizeof(unsigned char*));
if(bitp1==0||bitp2==0)
error("! no memory");

/*:16*//*22:*/
#line 376 "./hcoset.w"

for(int i= 0;i<FACT4;i++){
permtobit[i]= i;
bittoperm[i]= i;
}

/*:22*//*31:*/
#line 576 "./hcoset.w"

const int F2= 1+TWISTS;
const int R2= 1+2*TWISTS;
const int B2= 1+4*TWISTS;
const int L2= 1+5*TWISTS;
permcube pc;
for(int i= 0;i<FACT4/2;i++){
permtobit[2*i]= i;
pc.emp= 2*i;
pc.move(F2);
permtobit[pc.emp]= 12+i;
}
for(int i= 0;i<FACT4;i++)
bittoperm[permtobit[i]]= i;

/*:31*//*33:*/
#line 601 "./hcoset.w"

const int mvs[]= {R2,B2,L2};
for(int mvi= 0;mvi<SQMOVES;mvi++)
for(int p= 0;p<FACT4;p++){
pc.emp= p;
pc.move(mvs[mvi]);
rearrange[p&1][mvi][1<<(permtobit[p]%12)]= 
1<<(permtobit[pc.emp]%12);
}
for(int p= 0;p<2;p++)
for(int mvi= 0;mvi<SQMOVES;mvi++)
for(int i= 1;i<(1<<12);i++){
int lowb= i&-i;
rearrange[p][mvi][i]= 
rearrange[p][mvi][lowb]|rearrange[p][mvi][i-lowb];
}

/*:33*//*34:*/
#line 623 "./hcoset.w"

for(int i= 0;i<(1<<12);i++)
if(rearrange[0][1][i]!=rearrange[1][1][i])
error("! mismatch in rearrange");

/*:34*//*36:*/
#line 649 "./hcoset.w"

int ind= 0;
for(int e8_4= 0;e8_4<C8_4;e8_4++)
for(int epp1= 0;epp1<FACT4;epp1++)
for(int epp2= 0;epp2<FACT4;epp2+= 2,ind++){
int mvi= 0;
for(int mv= 0;mv<NMOVES;mv++){
if(!kocsymm::in_Kociemba_group(mv))
continue;
unpack_edgecoord(pc,e8_4,epp1,epp2);
pc.move(mv);
eperm_map[ind][mvi]= 
((permcube::c12_8[pc.et]*FACT4+pc.etp)*FACT4+pc.ebp)/2*3;
mvi++;
}
}

/*:36*//*40:*/
#line 787 "./hcoset.w"

#ifdef LEVELCOUNTS
for(int i= 1;i<(1<<12);i++)
bc[i]= 1+bc[i&(i-1)];
#endif

/*:40*//*47:*/
#line 922 "./hcoset.w"

unsigned char used[70];
memset(initorder,255,sizeof(initorder));
memset(used,255,sizeof(used));
int at= 0;
for(int i= 0;i<70;i++)
if(used[i]==255){
permcube pc2;
used[i]= 0;
initorder[at++]= i;
unpack_edgecoord(pc,i,0,0);
for(int j= 0;j<15;j++){
pc2= pc;
pc2.move(moveseq16[j]);
int e84= permcube::c12_8[pc.et];
if(used[e84]==255){
used[e84]= 0;
initorder[at++]= e84;
}
}
}
if(at!=70)
error("! bad setup in setorder()");

/*:47*//*50:*/
#line 970 "./hcoset.w"

pthread_mutex_init(&mutex,NULL);

/*:50*//*84:*/
#line 1839 "./hcoset.w"

#ifdef LEVELCOUNTS
setupk3map();
#endif

/*:84*/
#line 36 "./hcoset.w"

/*81:*/
#line 1759 "./hcoset.w"

if(argc> 1){
docoset(0,argv[1]);
}else{
if(coset_count==0)
coset_count= TOTALCOSETS;
genseqs(first_coset,first_coset+coset_count);
}

/*:81*/
#line 37 "./hcoset.w"

/*54:*/
#line 1027 "./hcoset.w"

if(singfile){
cout<<"Wrote "<<singcount<<" total positions to singfile"<<endl;
fclose(singfile);
}

/*:54*//*86:*/
#line 1895 "./hcoset.w"

#ifdef LEVELCOUNTS
for(unsigned int i= 0;i<sizeof(sum_ulev)/sizeof(sum_ulev[0]);i++)
if(sum_ulev[i])
cout<<"Level "<<i<<" count "<<sum_ulev[i]<<endl;
#endif/*:86*/
#line 38 "./hcoset.w"

phase1prune::check_integrity();
cout<<"Completed in "<<(walltime()-progstart)<<endl;
}

/*:1*/
