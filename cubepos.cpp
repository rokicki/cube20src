/*11:*/
#line 251 "./cubepos.w"

#include <iostream> 
#include "cubepos.h"
#include <math.h> 
/*13:*/
#line 279 "./cubepos.w"

const cubepos identity_cube(0,0,0);
unsigned char cubepos::corner_ori_inc[CUBIES],
cubepos::corner_ori_dec[CUBIES],
cubepos::corner_ori_neg_strip[CUBIES],cubepos::mod24[2*CUBIES];

/*:13*//*19:*/
#line 369 "./cubepos.w"

char cubepos::faces[FACES]= {'U','F','R','D','B','L'};

/*:19*//*22:*/
#line 413 "./cubepos.w"

unsigned char cubepos::edge_trans[NMOVES][CUBIES],
cubepos::corner_trans[NMOVES][CUBIES];

/*:22*//*24:*/
#line 451 "./cubepos.w"

static const unsigned char edge_twist_perm[FACES][4]= {
{0,2,3,1},
{3,7,11,6},
{2,5,10,7},
{9,11,10,8},
{0,4,8,5},
{1,6,9,4}
};

/*:24*//*25:*/
#line 470 "./cubepos.w"

static const unsigned char corner_twist_perm[FACES][4]= {
{0,1,3,2},
{2,3,7,6},
{3,1,5,7},
{4,6,7,5},
{1,0,4,5},
{0,2,6,4}
};

/*:25*//*26:*/
#line 515 "./cubepos.w"

static const unsigned char edge_change[FACES]= {0,0,1,0,0,1};

/*:26*//*27:*/
#line 546 "./cubepos.w"

static const unsigned char corner_change[FACES][4]= {
{0,0,0,0},
{1,2,1,2},
{1,2,1,2},
{0,0,0,0},
{1,2,1,2},
{1,2,1,2},
};

/*:27*//*33:*/
#line 643 "./cubepos.w"

unsigned char cubepos::inv_move[NMOVES];

/*:33*//*44:*/
#line 877 "./cubepos.w"

static char static_buf[200];

/*:44*//*47:*/
#line 959 "./cubepos.w"

static const char*sing_solved= 
"UF UR UB UL DF DR DB DL FR FL BR BL UFR URB UBL ULF DRF DFL DLB DBR";

/*:47*//*48:*/
#line 972 "./cubepos.w"

static const char*const smedges[]= {
"UB","BU","UL","LU","UR","RU","UF","FU",
"LB","BL","RB","BR","LF","FL","RF","FR",
"DB","BD","DL","LD","DR","RD","DF","FD",
};
static const char*const smcorners[]= {
"UBL","URB","ULF","UFR","DLB","DBR","DFL","DRF",
"LUB","BUR","FUL","RUF","BDL","RDB","LDF","FDR",
"BLU","RBU","LFU","FRU","LBD","BRD","FLD","RFD",
"ULB","UBR","UFL","URF","DBL","DRB","DLF","DFR",
"LBU","BRU","FLU","RFU","BLD","RBD","LFD","FRD",
"BUL","RUB","LUF","FUR","LDB","BDR","FDL","RDF",
};

/*:48*//*49:*/
#line 995 "./cubepos.w"

const int INVALID= 99;
static unsigned char lookup_edge_cubie[FACES*FACES];
static unsigned char lookup_corner_cubie[FACES*FACES*FACES];
static unsigned char sm_corner_order[8];
static unsigned char sm_edge_order[12];
static unsigned char sm_edge_flipped[12];

/*:49*//*56:*/
#line 1187 "./cubepos.w"

unsigned char cubepos::face_map[M][FACES],cubepos::move_map[M][NMOVES];
unsigned char cubepos::mm[M][M],cubepos::invm[M];
unsigned char cubepos::rot_edge[M][CUBIES],cubepos::rot_corner[M][CUBIES];

/*:56*//*57:*/
#line 1202 "./cubepos.w"

static const char*const axis_permute_map[]= 
{"UFR","URF","FRU","FUR","RUF","RFU"};
static const char*const axis_negate_map[]= 
{"UFR","UFL","UBL","UBR","DBR","DBL","DFL","DFR"};

/*:57*//*72:*/
#line 1458 "./cubepos.w"

unsigned char cubepos::canon_seq[CANONSEQSTATES][NMOVES];
int cubepos::canon_seq_mask[CANONSEQSTATES];
int cubepos::canon_seq_mask_ext[CANONSEQSTATES];

/*:72*/
#line 255 "./cubepos.w"

/*36:*/
#line 707 "./cubepos.w"

void cubepos::invert_into(cubepos&dst)const{
for(int i= 0;i<8;i++){
int cval= c[i];
dst.c[corner_perm(cval)]= corner_ori_sub(i,cval);
}
for(int i= 0;i<12;i++){
int cval= e[i];
dst.e[edge_perm(cval)]= edge_val(i,edge_ori(cval));
}
}

/*:36*//*50:*/
#line 1007 "./cubepos.w"

static int parse_cubie(const char*&p){
cubepos::skip_whitespace(p);
int v= 1;
int f= 0;
while((f= cubepos::parse_face(p))>=0){
v= v*6+f;
if(v>=2*6*6*6)
return-1;
}
return v;
}
static int parse_edge(const char*&p){
int c= parse_cubie(p);
if(c<6*6||c>=2*6*6)
return-1;
c= lookup_edge_cubie[c-6*6];
if(c==INVALID)
return-1;
return c;
}
static int parse_corner(const char*&p){
int c= parse_cubie(p);
if(c<6*6*6||c>=2*6*6*6)
return-1;
c= lookup_corner_cubie[c-6*6*6];
if(c==INVALID||c>=CUBIES)
return-1;
return c;
}

/*:50*//*53:*/
#line 1074 "./cubepos.w"

const char*cubepos::parse_Singmaster(const char*p){
if(strncmp(p,"SING ",5)==0)
p+= 5;
int m= 0;
for(int i= 0;i<12;i++){
int c= parse_edge(p)^sm_edge_flipped[i];
if(c<0)
return"No such edge";
e[edge_perm(c)]= edge_val(sm_edge_order[i],edge_ori(c));
m|= 1<<i;
}
for(int i= 0;i<8;i++){
int cval= parse_corner(p);
if(cval<0)
return"No such corner";
c[corner_perm(cval)]= corner_ori_sub(sm_corner_order[i],cval);
m|= 1<<(i+12);
}
skip_whitespace(p);
if(*p)
return"Extra stuff after Singmaster representation";
if(m!=((1<<20)-1))
return"Missing at least one cubie";
return 0;
}

/*:53*//*58:*/
#line 1212 "./cubepos.w"

static void parse_corner_to_facemap(const char*p,unsigned char*a){
for(int i= 0;i<3;i++){
int f= cubepos::parse_face(p[i]);
a[i]= f;
a[i+3]= (f+3)%FACES;
}
}

/*:58*//*59:*/
#line 1224 "./cubepos.w"

static void face_map_multiply(unsigned char*a,unsigned char*b,
unsigned char*c){
for(int i= 0;i<6;i++)
c[i]= b[a[i]];
}

/*:59*/
#line 256 "./cubepos.w"

void cubepos::init(){
static int initialized= 0;
if(initialized)
return;
initialized= 1;
/*14:*/
#line 287 "./cubepos.w"

for(int i= 0;i<CUBIES;i++){
int perm= corner_perm(i);
int ori= corner_ori(i);
corner_ori_inc[i]= corner_val(perm,(ori+1)%3);
corner_ori_dec[i]= corner_val(perm,(ori+2)%3);
corner_ori_neg_strip[i]= corner_val(0,(3-ori)%3);
mod24[i]= mod24[i+CUBIES]= i;
}

/*:14*//*28:*/
#line 563 "./cubepos.w"

for(int m= 0;m<NMOVES;m++)
for(int c= 0;c<CUBIES;c++){
edge_trans[m][c]= c;
corner_trans[m][c]= c;
}

/*:28*//*29:*/
#line 578 "./cubepos.w"

for(int f= 0;f<FACES;f++)
for(int t= 0;t<3;t++){
int m= f*TWISTS+t;
int isquarter= (t==0||t==2);
int perminc= t+1;
if(m<0)
continue;
for(int i= 0;i<4;i++){
int ii= (i+perminc)%4;
for(int o= 0;o<2;o++){
int oo= o;
if(isquarter)
oo^= edge_change[f];
edge_trans[m][edge_val(edge_twist_perm[f][i],o)]= 
edge_val(edge_twist_perm[f][ii],oo);
}
for(int o= 0;o<3;o++){
int oo= o;
if(isquarter)
oo= (corner_change[f][i]+oo)%3;
corner_trans[m][corner_val(corner_twist_perm[f][i],o)]= 
corner_val(corner_twist_perm[f][ii],oo);
}
}
}

/*:29*//*34:*/
#line 649 "./cubepos.w"

for(int i= 0;i<NMOVES;i++)
inv_move[i]= TWISTS*(i/TWISTS)+(NMOVES-i-1)%TWISTS;

/*:34*//*51:*/
#line 1040 "./cubepos.w"

memset(lookup_edge_cubie,INVALID,sizeof(lookup_edge_cubie));
memset(lookup_corner_cubie,INVALID,sizeof(lookup_corner_cubie));
for(int i= 0;i<CUBIES;i++){
const char*tmp= 0;
lookup_corner_cubie[parse_cubie(tmp= smcorners[i])-6*6*6]= i;
lookup_corner_cubie[parse_cubie(tmp= smcorners[CUBIES+i])-6*6*6]= CUBIES+i;
lookup_edge_cubie[parse_cubie(tmp= smedges[i])-6*6]= i;
}
const char*p= sing_solved;
for(int i= 0;i<12;i++){
int cv= parse_edge(p);
sm_edge_order[i]= edge_perm(cv);
sm_edge_flipped[i]= edge_ori(cv);
}
for(int i= 0;i<8;i++)
sm_corner_order[i]= corner_perm(parse_corner(p));

/*:51*//*60:*/
#line 1238 "./cubepos.w"

unsigned char face_to_m[FACES*FACES*FACES];
for(int i= 0;i<6;i++)
parse_corner_to_facemap(axis_permute_map[i],face_map[8*i]);
for(int i= 0;i<8;i++)
parse_corner_to_facemap(axis_negate_map[i],face_map[i]);
for(int i= 1;i<6;i++)
for(int j= 1;j<8;j++)
face_map_multiply(face_map[8*i],face_map[j],face_map[8*i+j]);

/*:60*//*61:*/
#line 1251 "./cubepos.w"

for(int i= 0;i<M;i++){
int v= face_map[i][0]*36+face_map[i][1]*6+face_map[i][2];
face_to_m[v]= i;
}
unsigned char tfaces[6];
for(int i= 0;i<M;i++)
for(int j= 0;j<M;j++){
face_map_multiply(face_map[i],face_map[j],tfaces);
int v= tfaces[0]*36+tfaces[1]*6+tfaces[2];
mm[i][j]= face_to_m[v];
if(mm[i][j]==0)
invm[i]= j;
}
for(int m= 0;m<M;m++){
int is_neg= (m^(m>>3))&1;
for(int f= 0;f<6;f++){
for(int t= 0;t<TWISTS;t++){
if(is_neg)
move_map[m][f*TWISTS+t]= face_map[m][f]*TWISTS+TWISTS-1-t;
else
move_map[m][f*TWISTS+t]= face_map[m][f]*TWISTS+t;
}
}
}

/*:61*//*62:*/
#line 1281 "./cubepos.w"

for(int m= 0;m<M;m++)
for(int c= 0;c<CUBIES;c++){
int v= 0;
for(int i= 0;i<2;i++)
v= 6*v+face_map[m][parse_face(smedges[c][i])];
rot_edge[m][c]= lookup_edge_cubie[v];
v= 0;
for(int i= 0;i<3;i++)
v= 6*v+face_map[m][parse_face(smcorners[c][i])];
rot_corner[m][c]= mod24[lookup_corner_cubie[v]];
}

/*:62*//*73:*/
#line 1467 "./cubepos.w"

for(int s= 0;s<CANONSEQSTATES;s++){
int prevface= (s-1)%FACES;
canon_seq_mask[s]= (1<<NMOVES)-1;
for(int mv= 0;mv<NMOVES;mv++){
int f= mv/TWISTS;
int isplus= 0;
if(s!=0&&(prevface==f||prevface==f+3))
{
canon_seq[s][mv]= INVALID;
canon_seq_mask[s]&= ~(1<<mv);
}else{
canon_seq[s][mv]= f+1+FACES*isplus;
}
}
canon_seq_mask_ext[s]= canon_seq_mask[s];
}

/*:73*/
#line 262 "./cubepos.w"

}

/*:11*//*17:*/
#line 328 "./cubepos.w"

cubepos::cubepos(int,int,int){
for(int i= 0;i<8;i++)
c[i]= corner_val(i,0);
for(int i= 0;i<12;i++)
e[i]= edge_val(i,0);
init();
}

/*:17*//*23:*/
#line 427 "./cubepos.w"

void cubepos::move(int mov){
const unsigned char*p= corner_trans[mov];
c[0]= p[c[0]];c[1]= p[c[1]];c[2]= p[c[2]];
c[3]= p[c[3]];c[4]= p[c[4]];c[5]= p[c[5]];
c[6]= p[c[6]];c[7]= p[c[7]];
p= edge_trans[mov];
e[0]= p[e[0]];e[1]= p[e[1]];e[2]= p[e[2]];
e[3]= p[e[3]];e[4]= p[e[4]];e[5]= p[e[5]];
e[6]= p[e[6]];e[7]= p[e[7]];e[8]= p[e[8]];
e[9]= p[e[9]];e[10]= p[e[10]];e[11]= p[e[11]];
}

/*:23*//*35:*/
#line 658 "./cubepos.w"

moveseq cubepos::invert_sequence(const moveseq&seq){
unsigned int len= seq.size();
moveseq r(len);
for(unsigned int i= 0;i<len;i++)
r[len-i-1]= invert_move(seq[i]);
return r;
}

/*:35*//*38:*/
#line 761 "./cubepos.w"

#define ROT2(cc,a,b) {unsigned char t= cc[a];cc[a]= cc[b];cc[b]= t;}
#define ROT4(cc,a,b,c,d) {unsigned char t= cc[d];cc[d]= cc[c];cc[c]= cc[b]; cc[b]= cc[a];cc[a]= t;}
#define ROT22(cc,a,b,c,d) ROT2(cc,a,c) ROT2(cc,b,d)

/*:38*//*39:*/
#line 773 "./cubepos.w"

#define EDGE4FLIP(a,b,c,d) {unsigned char t= e[d];e[d]= edge_flip(e[c]);\
         e[c]= edge_flip(e[b]); e[b]= edge_flip(e[a]);e[a]= edge_flip(t);}
#define CORNER4FLIP(a,b,cc,d) {unsigned char t= c[d];c[d]= corner_ori_inc[c[cc]];\
 c[cc]= corner_ori_dec[c[b]];c[b]= corner_ori_inc[c[a]];c[a]= corner_ori_dec[t];}

/*:39*//*40:*/
#line 788 "./cubepos.w"

void cubepos::movepc(int mov){
switch(mov){
case 0:ROT4(e,0,2,3,1);ROT4(c,0,1,3,2);break;
case 1:ROT22(e,0,2,3,1);ROT22(c,0,1,3,2);break;
case 2:ROT4(e,1,3,2,0);ROT4(c,2,3,1,0);break;
case 3:ROT4(e,3,7,11,6);CORNER4FLIP(3,7,6,2);break;
case 4:ROT22(e,3,7,11,6);ROT22(c,2,3,7,6);break;
case 5:ROT4(e,6,11,7,3);CORNER4FLIP(3,2,6,7);break;
case 6:EDGE4FLIP(2,5,10,7);CORNER4FLIP(1,5,7,3);break;
case 7:ROT22(e,2,5,10,7);ROT22(c,3,1,5,7);break;
case 8:EDGE4FLIP(7,10,5,2);CORNER4FLIP(1,3,7,5);break;
case 9:ROT4(e,9,11,10,8);ROT4(c,4,6,7,5);break;
case 10:ROT22(e,9,11,10,8);ROT22(c,4,6,7,5);break;
case 11:ROT4(e,8,10,11,9);ROT4(c,5,7,6,4);break;
case 12:ROT4(e,0,4,8,5);CORNER4FLIP(0,4,5,1);break;
case 13:ROT22(e,0,4,8,5);ROT22(c,1,0,4,5);break;
case 14:ROT4(e,5,8,4,0);CORNER4FLIP(0,1,5,4);break;
case 15:EDGE4FLIP(1,6,9,4);CORNER4FLIP(2,6,4,0);break;
case 16:ROT22(e,1,6,9,4);ROT22(c,0,2,6,4);break;
case 17:EDGE4FLIP(4,9,6,1);CORNER4FLIP(2,0,4,6);break;
}
}

/*:40*//*42:*/
#line 845 "./cubepos.w"

void cubepos::mul(const cubepos&a,const cubepos&b,cubepos&r){
for(int i= 0;i<8;i++){
int cc= a.c[i];
r.c[i]= corner_ori_add(b.c[corner_perm(cc)],cc);
}
for(int i= 0;i<12;i++){
int cc= a.e[i];
r.e[i]= edge_ori_add(b.e[edge_perm(cc)],cc);
}
}

/*:42*//*45:*/
#line 882 "./cubepos.w"

void cubepos::skip_whitespace(const char*&p){
while(*p&&*p<=' ')
p++;
}
int cubepos::parse_face(const char*&p){
int f= parse_face(*p);
if(f>=0)
p++;
return f;
}
int cubepos::parse_face(char f){
switch(f){
case'u':case'U':return 0;
case'f':case'F':return 1;
case'r':case'R':return 2;
case'd':case'D':return 3;
case'b':case'B':return 4;
case'l':case'L':return 5;
default:
return-1;
}
}
int cubepos::parse_move(const char*&p){
skip_whitespace(p);
const char*q= p;
int f= parse_face(q);
if(f<0)
return-1;
int t= 0;
switch(*q){
case'1':case'+':t= 0;break;
case'2':t= 1;break;
case'3':case'\'':case'-':t= TWISTS-1;break;
default:
return-1;
}
p= q+1;
return f*TWISTS+t;
}

/*:45*//*46:*/
#line 925 "./cubepos.w"

void cubepos::append_move(char*&p,int mv){
append_face(p,mv/TWISTS);
*p++= "123"[mv%TWISTS];
*p= 0;
}
moveseq cubepos::parse_moveseq(const char*&p){
moveseq r;
int mv;
while((mv= parse_move(p))>=0)
r.push_back(mv);
return r;
}
void cubepos::append_moveseq(char*&p,const moveseq&seq){
*p= 0;
for(unsigned int i= 0;i<seq.size();i++)
append_move(p,seq[i]);
}
char*cubepos::moveseq_string(const moveseq&seq){
if(seq.size()> 65)
error("! can't print a move sequence that long");
char*p= static_buf;
append_moveseq(p,seq);
return static_buf;
}

/*:46*//*54:*/
#line 1106 "./cubepos.w"

char*cubepos::Singmaster_string()const{
cubepos cp;
invert_into(cp);
char*p= static_buf;
for(int i= 0;i<12;i++){
if(i!=0)
*p++= ' ';
int j= sm_edge_order[i];
const char*q= smedges[cp.e[j]^sm_edge_flipped[i]];
*p++= *q++;
*p++= *q++;
}
for(int i= 0;i<8;i++){
*p++= ' ';
int j= sm_corner_order[i];
const char*q= smcorners[cp.c[j]];
*p++= *q++;
*p++= *q++;
*p++= *q++;
}
*p= 0;
return static_buf;
}

/*:54*//*65:*/
#line 1316 "./cubepos.w"

void cubepos::remap_into(int m,cubepos&dst)const{
int mprime= invm[m];
for(int i= 0;i<8;i++){
int c1= rot_corner[mprime][i];
int c2= corner_ori_add(c[corner_perm(c1)],c1);
dst.c[i]= rot_corner[m][c2];
}
for(int i= 0;i<12;i++){
int c1= rot_edge[mprime][i*2];
int c2= edge_ori_add(e[edge_perm(c1)],c1);
dst.e[i]= rot_edge[m][c2];
}
}

/*:65*//*66:*/
#line 1335 "./cubepos.w"

void cubepos::canon_into48_aux(cubepos&dst)const{
for(int m= 1;m<M;m++){
int mprime= invm[m];
int isless= 0;
for(int i= 0;i<8;i++){
int c1= rot_corner[mprime][i];
int c2= corner_ori_add(c[corner_perm(c1)],c1);
int t= rot_corner[m][c2];
if(isless||t<dst.c[i]){
dst.c[i]= t;
isless= 1;
}else if(t> dst.c[i])
goto nextm;
}
for(int i= 0;i<12;i++){
int c1= rot_edge[mprime][i*2];
int c2= edge_ori_add(e[edge_perm(c1)],c1);
int t= rot_edge[m][c2];
if(isless||t<dst.e[i]){
dst.e[i]= t;
isless= 1;
}else if(t> dst.e[i])
goto nextm;
}
nextm:;
}
}
void cubepos::canon_into48(cubepos&dst)const{
dst= *this;
canon_into48_aux(dst);
}

/*:66*//*68:*/
#line 1383 "./cubepos.w"

void cubepos::randomize(){
int parity= 0;
for(int i= 0;i<7;i++){
int j= i+(int)((8-i)*myrand());
if(i!=j){
swap(c[i],c[j]);
parity++;
}
}
for(int i= 0;i<11;i++){
int j= i+(int)((12-i)*myrand());
if(i!=j){
swap(e[i],e[j]);
parity++;
}
}
if(parity&1)
swap(e[10],e[11]);
int s= 24;
for(int i= 0;i<7;i++){
int a= (int)(3*myrand());
s-= a;
c[i]= corner_val(corner_perm(c[i]),a);
}
c[7]= corner_val(corner_perm(c[7]),s%3);
s= 0;
for(int i= 0;i<11;i++){
int a= (int)(2*myrand());
e[i]= edge_ori_add(e[i],a);
s^= a;
}
e[11]^= s;
}

/*:68*//*69:*/
#line 1422 "./cubepos.w"

void cubepos::canon_into96(cubepos&dst)const{
cubepos cpi;
invert_into(cpi);
if(*this<cpi){
dst= *this;
}else{
dst= cpi;
}
canon_into48_aux(dst);
cpi.canon_into48_aux(dst);
}

/*:69*//*76:*/
#line 1506 "./cubepos.w"

void error(const char*s){
cerr<<s<<endl;
if(*s=='!')
exit(10);
}
static double start;
double walltime(){
struct timeval tv;
gettimeofday(&tv,0);
return tv.tv_sec+0.000001*tv.tv_usec;
}
double duration(){
double now= walltime();
double r= now-start;
start= now;
return r;
}


/*:76*/
