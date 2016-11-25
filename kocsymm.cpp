/*7:*/
#line 142 "./kocsymm.w"

#include "kocsymm.h"
#include <iostream> 
using namespace std;
/*8:*/
#line 160 "./kocsymm.w"

lookup_type kocsymm::cornermove[CORNERSYMM][NMOVES];
lookup_type kocsymm::edgeomove[EDGEOSYMM][NMOVES];
lookup_type kocsymm::edgepmove[EDGEPERM][NMOVES];

/*:8*//*12:*/
#line 214 "./kocsymm.w"

lookup_type kocsymm::epsymm_compress[1<<12];
lookup_type kocsymm::epsymm_expand[EDGEOSYMM];

/*:12*//*21:*/
#line 395 "./kocsymm.w"

lookup_type kocsymm::cornersymm_expand[CORNERRSYMM];
corner_mapinfo kocsymm::cornersymm[CORNERSYMM];
lookup_type kocsymm::edgeomap[EDGEOSYMM][KOCSYMM];
lookup_type kocsymm::edgepmap[EDGEPERM][KOCSYMM];
lookup_type kocsymm::edgepxor[EDGEPERM][2];

/*:21*//*32:*/
#line 588 "./kocsymm.w"

unsigned char permcube::s4inv[FACT4];
unsigned char permcube::s4mul[FACT4][FACT4];
unsigned char permcube::s4compress[256];
unsigned char permcube::s4expand[FACT4];

/*:32*//*37:*/
#line 675 "./kocsymm.w"

unsigned char permcube::c8_4_compact[256];
unsigned char permcube::c8_4_expand[C8_4];
unsigned char permcube::c8_4_parity[C8_4];

/*:37*//*40:*/
#line 713 "./kocsymm.w"

unsigned char permcube::c12_8[EDGEPERM];
lookup_type permcube::c8_12[C8_4];

/*:40*//*44:*/
#line 778 "./kocsymm.w"

unsigned short permcube::eperm_move[EDGEPERM][NMOVES];
int permcube::cperm_move[C8_4][NMOVES];

/*:44*/
#line 146 "./kocsymm.w"

/*13:*/
#line 222 "./kocsymm.w"

static int bc(int v){
int r= 0;
while(v){
v&= v-1;
r++;
}
return r;
}

/*:13*//*35:*/
#line 655 "./kocsymm.w"

int muls4(int a,int b){
int r= 3&(b>>(2*(a&3)));
r+= (3&(b>>(2*((a>>2)&3))))<<2;
r+= (3&(b>>(2*((a>>4)&3))))<<4;
r+= (3&(b>>(2*((a>>6)&3))))<<6;
return r;
}

/*:35*/
#line 147 "./kocsymm.w"

/*15:*/
#line 254 "./kocsymm.w"

kocsymm::kocsymm(const cubepos&cp){
int c= 0,eo= 0,ep= 0;
for(int i= 6;i>=0;i--)
c= 3*c+cubepos::corner_ori(cp.c[i]);
for(int i= 10;i>=0;i--){
eo= 2*eo+cubepos::edge_ori(cp.e[i]);
ep= 2*ep+(cp.e[i]&8);
}
csymm= c;
eosymm= eo;
epsymm= epsymm_compress[ep>>3];
}

/*:15*//*16:*/
#line 274 "./kocsymm.w"

void kocsymm::set_coset(cubepos&cp){
int c= csymm,eo= eosymm,ep= epsymm_expand[epsymm];
int s= 0;
for(int i= 0;i<7;i++){
int ori= c%3;
cp.c[i]= cubepos::corner_val(i,ori);
s+= ori;
c= c/3;
}
cp.c[7]= cubepos::corner_val(7,(8*3-s)%3);
s= 0;
int nextmid= 4;
int nextud= 0;
for(int i= 0;i<12;i++){
if(i==11)
eo= s;
int ori= eo&1;
if(ep&1)
cp.e[i]= cubepos::edge_val(nextmid++,ori);
else{
cp.e[i]= cubepos::edge_val(nextud++,ori);
if(nextud==4)
nextud= 8;
}
s^= ori;
eo>>= 1;
ep>>= 1;
}
}

/*:16*//*25:*/
#line 474 "./kocsymm.w"

void kocsymm::canon_into(kocsymm&kc)const{
corner_mapinfo&cm= cornersymm[csymm];
kc.csymm= cornersymm_expand[cm.csymm];
kc.eosymm= edgeomap[edgepxor[epsymm][cm.minmap>>3]^eosymm][cm.minmap];
kc.epsymm= edgepmap[epsymm][cm.minmap];
for(int m= cm.minmap+1;cm.minbits>>m;m++)
if((cm.minbits>>m)&1){
int neo= edgeomap[edgepxor[epsymm][m>>3]^eosymm][m];
if(neo> kc.eosymm)
continue;
int nep= edgepmap[epsymm][m];
if(neo<kc.eosymm||nep<kc.epsymm){
kc.eosymm= neo;
kc.epsymm= nep;
}
}
}

/*:25*//*27:*/
#line 500 "./kocsymm.w"

int kocsymm::calc_symm()const{
int r= 1;
corner_mapinfo&cm= cornersymm[csymm];
int teosymm= edgeomap[edgepxor[epsymm][cm.minmap>>3]^eosymm][cm.minmap];
int tepsymm= edgepmap[epsymm][cm.minmap];
for(int m= cm.minmap+1;cm.minbits>>m;m++)
if(((cm.minbits>>m)&1)&&
edgeomap[edgepxor[epsymm][m>>3]^eosymm][m]==teosymm&&
edgepmap[epsymm][m]==tepsymm)
r++;
return r;
}

/*:27*//*33:*/
#line 598 "./kocsymm.w"

void permcube::init(){
/*34:*/
#line 632 "./kocsymm.w"

int cc= 0;
for(int a= 0;a<4;a++)
for(int b= 0;b<4;b++)if(a!=b)
for(int c= 0;c<4;c++)if(a!=c&&b!=c){
int d= 0+1+2+3-a-b-c;
int coor= cc^((cc>>1)&1);
int expanded= (1<<(2*b))+(2<<(2*c))+(3<<(2*d));
s4compress[expanded]= coor;
s4expand[coor]= expanded;
cc++;
}
for(int i= 0;i<FACT4;i++)
for(int j= 0;j<FACT4;j++){
int k= s4compress[muls4(s4expand[i],s4expand[j])];
s4mul[j][i]= k;
if(k==0)
s4inv[i]= j;
}

/*:34*//*38:*/
#line 684 "./kocsymm.w"

int c= 0;
for(int i= 0;i<256;i++)
if(bc(i)==4){
int parity= 0;
for(int j= 0;j<8;j++)
if(1&(i>>j))
for(int k= 0;k<j;k++)
if(0==(1&(i>>k)))
parity++;
c8_4_parity[c]= parity&1;
c8_4_compact[i]= c;
c8_4_expand[c]= i;
c++;
}

/*:38*//*41:*/
#line 720 "./kocsymm.w"

for(int i= 0;i<EDGEPERM;i++){
int expbits= kocsymm::epsymm_expand[i];
if(expbits&0x0f0)
c12_8[i]= 255;
else{
int ii= c8_4_compact[(expbits>>4)+(expbits&15)];
c12_8[i]= ii;
c8_12[ii]= i;
}
}

/*:41*//*53:*/
#line 962 "./kocsymm.w"

cubepos cp,cp2;
for(int i= 0;i<EDGEPERM;i++){
permcube pc;
pc.em= i;
int remaining_edges= 0xfff-kocsymm::epsymm_expand[i];
int mask= 0;
int bitsseen= 0;
while(bitsseen<4){
if(remaining_edges&(mask+1))
bitsseen++;
mask= 2*mask+1;
}
pc.et= kocsymm::epsymm_compress[remaining_edges&mask];
pc.eb= kocsymm::epsymm_compress[remaining_edges&~mask];
pc.set_perm(cp);
for(int mv= 0;mv<NMOVES;mv++){
cp2= cp;
cp2.movepc(mv);
permcube pc2(cp2);
eperm_move[i][mv]= (pc2.em<<5)+pc2.emp;
}
}

/*:53*//*54:*/
#line 991 "./kocsymm.w"

for(int i= 0;i<C8_4;i++){
permcube pc;
pc.c8_4= i;
pc.set_perm(cp);
for(int mv= 0;mv<NMOVES;mv++){
cp2= cp;
cp2.movepc(mv);
permcube pc2(cp2);
cperm_move[i][mv]= (pc2.c8_4<<10)+(pc2.ctp<<5)+pc2.cbp;
}
}

/*:54*/
#line 600 "./kocsymm.w"
;
}

/*:33*//*46:*/
#line 791 "./kocsymm.w"

void permcube::move(int mv){
#ifdef SAFETY_CHECKS
if((kocsymm::epsymm_expand[et]|kocsymm::epsymm_expand[em]|
kocsymm::epsymm_expand[eb])!=0xfff)
error("! bad pc in move");
#endif
int t= eperm_move[et][mv];
et= t>>5;
etp= s4mul[etp][t&31];
t= eperm_move[em][mv];
em= t>>5;
emp= s4mul[emp][t&31];
t= eperm_move[eb][mv];
eb= t>>5;
ebp= s4mul[ebp][t&31];
t= cperm_move[c8_4][mv];
c8_4= t>>10;
ctp= s4mul[ctp][(t>>5)&31];
cbp= s4mul[cbp][t&31];
}

/*:46*//*48:*/
#line 843 "./kocsymm.w"

void permcube::init_edge_from_cp(const cubepos&cp){
et= em= eb= 0;
etp= emp= ebp= 0;
for(int i= 11;i>=0;i--){
int perm= cubepos::edge_perm(cp.e[i]);
if(perm&4){
em|= 1<<i;
emp= 4*emp+(perm&3);
}else if(perm&8){
eb|= 1<<i;
ebp= 4*ebp+(perm&3);
}else{
et|= 1<<i;
etp= 4*etp+(perm&3);
}
}
et= kocsymm::epsymm_compress[et];
em= kocsymm::epsymm_compress[em];
eb= kocsymm::epsymm_compress[eb];
etp= s4compress[etp];
emp= s4compress[emp];
ebp= s4compress[ebp];
}

/*:48*//*49:*/
#line 870 "./kocsymm.w"

void permcube::init_corner_from_cp(const cubepos&cp){
c8_4= 0;
ctp= cbp= 0;
for(int i= 7;i>=0;i--){
int perm= cubepos::corner_perm(cp.c[i]);
if(perm&4){
cbp= 4*cbp+(perm&3);
}else{
c8_4|= 1<<i;
ctp= 4*ctp+(perm&3);
}
}
c8_4= c8_4_compact[c8_4];
ctp= s4compress[ctp];
cbp= s4compress[cbp];
}
permcube::permcube(const cubepos&cp){
init_edge_from_cp(cp);
init_corner_from_cp(cp);
}

/*:49*//*50:*/
#line 895 "./kocsymm.w"

void permcube::set_edge_perm(cubepos&cp)const{
int et_bits= kocsymm::epsymm_expand[et];
int em_bits= kocsymm::epsymm_expand[em];
int et_perm= s4expand[etp];
int em_perm= s4expand[emp];
int eb_perm= s4expand[ebp];
for(int i= 0;i<12;i++)
if((et_bits>>i)&1){
cp.e[i]= cubepos::edge_val((3&et_perm),
cubepos::edge_ori(cp.e[i]));
et_perm>>= 2;
}else if((em_bits>>i)&1){
cp.e[i]= cubepos::edge_val((3&em_perm)+4,
cubepos::edge_ori(cp.e[i]));
em_perm>>= 2;
}else{
cp.e[i]= cubepos::edge_val((3&eb_perm)+8,
cubepos::edge_ori(cp.e[i]));
eb_perm>>= 2;
}
}

/*:50*//*51:*/
#line 920 "./kocsymm.w"

void permcube::set_corner_perm(cubepos&cp)const{
int c8_4_bits= c8_4_expand[c8_4];
int ct_perm= s4expand[ctp];
int cb_perm= s4expand[cbp];
for(int i= 0;i<8;i++)
if((c8_4_bits>>i)&1){
cp.c[i]= cubepos::corner_val((3&ct_perm),
cubepos::corner_ori(cp.c[i]));
ct_perm>>= 2;
}else{
cp.c[i]= cubepos::corner_val((3&cb_perm)+4,
cubepos::corner_ori(cp.c[i]));
cb_perm>>= 2;
}
}
void permcube::set_perm(cubepos&cp)const{
set_edge_perm(cp);
set_corner_perm(cp);
}

/*:51*//*52:*/
#line 948 "./kocsymm.w"

permcube::permcube(){
c8_4= 0;
ctp= cbp= 0;
et= kocsymm::epsymm_compress[0xf];
em= 0;
eb= kocsymm::epsymm_compress[0xf00];
etp= emp= ebp= 0;
}

/*:52*/
#line 148 "./kocsymm.w"

void kocsymm::init(){
static int initialized= 0;
if(initialized)
return;
initialized= 1;
/*14:*/
#line 236 "./kocsymm.w"

int c= 0;
for(int i= 0;i<1<<12;i++)
if(bc(i)==4){
int rotval= ((i<<4)+(i>>8))&0xfff;
epsymm_compress[rotval]= c;
epsymm_compress[rotval&0x7ff]= c;
epsymm_expand[c]= rotval;
c++;
}

/*:14*//*17:*/
#line 311 "./kocsymm.w"

cubepos cp,cp2;
for(int i= 0;i<CORNERSYMM;i++){
kocsymm kc(i,i%EDGEOSYMM,i%EDGEPERM);
kc.set_coset(cp);
for(int mv= 0;mv<NMOVES;mv++){
cp2= cp;
cp2.movepc(mv);
kocsymm kc2(cp2);
cornermove[i][mv]= kc2.csymm;
if(i<EDGEOSYMM)
edgeomove[i][mv]= kc2.eosymm;
if(i<EDGEPERM)
edgepmove[i][mv]= kc2.epsymm;
}
}

/*:17*//*22:*/
#line 407 "./kocsymm.w"

c= 0;
for(int cs= 0;cs<CORNERSYMM;cs++){
int minval= cs;
int lowm= 0;
int lowbits= 1;
kocsymm kc(cs,0,0);
for(int m= 1;m<KOCSYMM;m++){
kc.set_coset(cp);
cp.remap_into(m,cp2);
kocsymm kc2(cp2);
if(kc2.csymm<minval){
minval= kc2.csymm;
lowbits= 1<<m;
lowm= m;
}else if(kc2.csymm==minval){
lowbits|= 1<<m;
}
}
if(minval==cs){
cornersymm_expand[c]= minval;
cornersymm[cs].csymm= c++;
}
cornersymm[cs].minbits= lowbits;
cornersymm[cs].minmap= lowm;
cornersymm[cs].csymm= cornersymm[minval].csymm;
}
if(c!=CORNERRSYMM)
error("! bad cornersym result");

/*:22*//*23:*/
#line 442 "./kocsymm.w"

for(int ep= 0;ep<EDGEPERM;ep++){
kocsymm kc(0,0,ep);
for(int m= 0;m<KOCSYMM;m++){
kc.set_coset(cp);
cp.remap_into(m,cp2);
kocsymm kc2(cp2);
edgepmap[ep][m]= kc2.epsymm;
if(m==8){
edgepxor[kc2.epsymm][0]= 0;
edgepxor[kc2.epsymm][1]= kc2.eosymm;
}
}
}
for(int eo= 0;eo<EDGEOSYMM;eo++){
kocsymm kc(0,eo,0);
for(int m= 0;m<KOCSYMM;m++){
kc.set_coset(cp);
cp.remap_into(m,cp2);
kocsymm kc2(cp2);
edgeomap[eo][m]= kc2.eosymm;
}
}

/*:23*/
#line 154 "./kocsymm.w"

permcube::init();
}

/*:7*/
