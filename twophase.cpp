/*1:*/
#line 13 "./twophase.w"

const char*BANNER= 
"This is twophase 1.0, (C) 2010 Tomas Rokicki.  All Rights Reserved.";
#include "phase1prune.h"
#include "phase2prune.h"
#include <pthread.h> 
#include <iostream> 
#include <map> 
using namespace std;
/*2:*/
#line 43 "./twophase.w"

int verbose= 1;
int numthreads= 1;
const int MAX_THREADS= 32;

/*:2*//*4:*/
#line 80 "./twophase.w"

int target_length= 0;
long long phase2limit= 0xffffffffffffffLL;
long long phase2total= 0LL;

/*:4*//*7:*/
#line 117 "./twophase.w"

int skipwrite= 0;

/*:7*//*15:*/
#line 202 "./twophase.w"

pthread_mutex_t mutex;

/*:15*//*19:*/
#line 249 "./twophase.w"

int axesmask= 63;

/*:19*//*28:*/
#line 438 "./twophase.w"

class solution{
public:
solution(const cubepos&cparg,int seqarg,long long p2parg,
moveseq&solarg){
cp= cparg;
seq= seqarg;
phase2probes= p2parg;
moves= solarg;
}
solution(){}
cubepos cp;
int seq;
long long phase2probes;
moveseq moves;
};
map<int,solution> queue;
int next_sequence= 1;
int missed_target= 0;
int solved= 0;

/*:28*/
#line 22 "./twophase.w"

/*10:*/
#line 140 "./twophase.w"

const int MAX_MOVES= 32;

/*:10*//*17:*/
#line 212 "./twophase.w"

void get_global_lock(){
pthread_mutex_lock(&mutex);
}
void release_global_lock(){
pthread_mutex_unlock(&mutex);
}


/*:17*//*22:*/
#line 307 "./twophase.w"

int sloweq(const cubepos&cp1,const cubepos&cp2){
cubepos cp3;
for(int m= 0;m<KOCSYMM;m++){
cp2.remap_into(m,cp3);
if(cp1==cp3)
return 1;
}
return 0;
}

/*:22*//*27:*/
#line 416 "./twophase.w"

void display(const cubepos&cp,int seq,long long phase2probes,moveseq sol){
phase2total+= phase2probes;
if(verbose||(int)sol.size()> target_length){
if((int)sol.size()> target_length)
cout<<"WARNING: missed target for "<<cp.Singmaster_string()
<<endl;
cout<<"Solution "<<seq<<" len "<<sol.size()<<" probes "
<<phase2probes<<endl;
cout<<cubepos::moveseq_string(sol)<<endl;
}
}

/*:27*//*29:*/
#line 461 "./twophase.w"

void report(const cubepos&cp,int seq,long long phase2probes,moveseq sol){
get_global_lock();
solved++;
if((int)sol.size()> target_length&&target_length)
missed_target++;
if(seq==next_sequence){
display(cp,seq,phase2probes,sol);
next_sequence++;
while(queue.find(next_sequence)!=queue.end()){
solution&s= queue[next_sequence];
display(s.cp,s.seq,s.phase2probes,s.moves);
queue.erase(next_sequence);
next_sequence++;
}
}else{
queue[seq]= solution(cp,seq,phase2probes,sol);
}
release_global_lock();
}

/*:29*//*32:*/
#line 514 "./twophase.w"

int getwork(cubepos&cp){
static int input_seq= 1;
const int BUFSIZE= 1000;
char buf[BUFSIZE+1];
get_global_lock();
if(fgets(buf,BUFSIZE,stdin)==0){
release_global_lock();
return-1;
}
if(cp.parse_Singmaster(buf)!=0){
cp= identity_cube;
const char*p= buf;
moveseq ms= cubepos::parse_moveseq(p);
if(*p)
error("! could not parse position");
for(unsigned int i= 0;i<ms.size();i++)
cp.move(ms[i]);
}
int r= input_seq++;
release_global_lock();
return r;
}

/*:32*/
#line 23 "./twophase.w"

/*11:*/
#line 155 "./twophase.w"

class twophasesolver{
public:
twophasesolver(){}
cubepos pos;
long long phase2probes;
int bestsol;
int finished;
int curm;
int solmap;
int seq;
/*12:*/
#line 174 "./twophase.w"

unsigned char moves[MAX_MOVES];
unsigned char bestmoves[MAX_MOVES];

/*:12*//*18:*/
#line 237 "./twophase.w"

kocsymm kc6[6],kccanon6[6];
cubepos cp6[6];
permcube pc6[6];
int mindepth[6];
char uniq[6];
int minmindepth;

/*:18*/
#line 166 "./twophase.w"

/*13:*/
#line 182 "./twophase.w"

void solve(int seqarg,cubepos&cp){
/*14:*/
#line 192 "./twophase.w"

pos= cp;
phase2probes= 0;
bestsol= MAX_MOVES;
finished= 0;
seq= seqarg;

/*:14*/
#line 184 "./twophase.w"

/*21:*/
#line 269 "./twophase.w"

minmindepth= MAX_MOVES;
cubepos cpi,cp2;
pos.invert_into(cpi);
int ind= 0;
for(int inv= 0;inv<2;inv++)
for(int mm= 0;mm<3;mm++,ind++){
int m= KOCSYMM*mm;
if(inv)
cpi.remap_into(m,cp2);
else
pos.remap_into(m,cp2);
cp6[ind]= cp2;
kc6[ind]= kocsymm(cp2);
pc6[ind]= permcube(cp2);
kc6[ind].canon_into(kccanon6[ind]);
mindepth[ind]= phase1prune::lookup(kc6[ind]);
if(mindepth[ind]<minmindepth)
minmindepth= mindepth[ind];
uniq[ind]= 1&(axesmask>>ind);
for(int i= 0;i<ind;i++)
if(uniq[i]&&kccanon6[ind]==kccanon6[i]&&
sloweq(cp6[ind],cp6[i])){
uniq[ind]= 0;
break;
}
if(verbose> 1){
get_global_lock();
cout<<"Axis "<<ind<<" depth "<<mindepth[ind]<<" uniq "
<<(int)uniq[ind]<<endl;
release_global_lock();
}
}

/*:21*/
#line 185 "./twophase.w"

/*23:*/
#line 320 "./twophase.w"

for(int d= minmindepth;d<bestsol&&!finished;d++){
for(curm= 0;curm<6;curm++)
if(uniq[curm]&&d<bestsol&&!finished&&d>=mindepth[curm]){
if(verbose> 1){
get_global_lock();
cout<<"Orientation "<<curm<<" at depth "<<d<<endl;
release_global_lock();
}
solvep1(kc6[curm],pc6[curm],d,0,ALLMOVEMASK,CANONSEQSTART);
}
}

/*:23*/
#line 186 "./twophase.w"

/*26:*/
#line 398 "./twophase.w"

moveseq sol;
int m= cubepos::invm[(solmap%3)*KOCSYMM];
for(int i= 0;i<bestsol;i++)
sol.push_back(cubepos::move_map[m][bestmoves[i]]);
if(solmap>=3)
sol= cubepos::invert_sequence(sol);
cubepos cpt;
for(unsigned int i= 0;i<sol.size();i++)
cpt.move(sol[i]);
if(cpt!=pos)
error("! move sequence doesn't work");
report(pos,seq,phase2probes,sol);

/*:26*/
#line 187 "./twophase.w"

}

/*:13*//*24:*/
#line 335 "./twophase.w"

void solvep1(const kocsymm&kc,const permcube&pc,int togo,int sofar,
int movemask,int canon){
if(togo==0){
if(kc==identity_kc)
solvep2(pc,sofar);
return;
}
if(finished)
return;
togo--;
kocsymm kc2;
permcube pc2;
int newmovemask;
while(!finished&&movemask){
int mv= ffs(movemask)-1;
movemask&= movemask-1;
kc2= kc;
kc2.move(mv);
int nd= phase1prune::lookup(kc2,togo,newmovemask);
if(nd<=togo&&(togo==nd||togo+nd>=5)){
pc2= pc;
pc2.move(mv);
moves[sofar]= mv;
int new_canon= cubepos::next_cs(canon,mv);
solvep1(kc2,pc2,togo,sofar+1,
newmovemask&cubepos::cs_mask(new_canon),
new_canon);
}
}
}

/*:24*//*25:*/
#line 369 "./twophase.w"

void solvep2(const permcube&pc,int sofar){
phase2probes++;
int d= phase2prune::lookup(pc);
if(d+sofar<bestsol){
moveseq ms= phase2prune::solve(pc,bestsol-sofar-1);
if((int)(ms.size())+sofar<bestsol&&
(ms.size()> 0||pc==identity_pc)){
bestsol= ms.size()+sofar;
for(unsigned int i= 0;i<ms.size();i++)
moves[sofar+i]= ms[i];
memcpy(bestmoves,moves,bestsol);
if(verbose> 1){
get_global_lock();
cout<<"New solution for "<<seq<<" at "<<bestsol<<endl;
release_global_lock();
}
solmap= curm;
if(bestsol<=target_length)
finished= 1;
}
}
if(phase2probes>=phase2limit&&bestsol<MAX_MOVES)
finished= 1;
}

/*:25*//*31:*/
#line 499 "./twophase.w"

void dowork(){
cubepos cp;
int seq;
while(1){
seq= getwork(cp);
if(seq<=0)
return;
solve(seq,cp);
}
}

/*:31*//*33:*/
#line 541 "./twophase.w"

static void*worker(void*s){
twophasesolver*solv= (twophasesolver*)s;
solv->dowork();
return 0;
}

/*:33*/
#line 167 "./twophase.w"

char pad[256];
}solvers[MAX_THREADS];

/*:11*/
#line 24 "./twophase.w"

int main(int argc,char*argv[]){
double progstart= walltime();
duration();
/*3:*/
#line 50 "./twophase.w"

while(argc> 1&&argv[1][0]=='-'){
argc--;
argv++;
switch(argv[0][1]){
case'v':verbose++;break;
case'q':verbose= 0;break;
case't':
if(sscanf(argv[1],"%d",&numthreads)!=1)
error("! bad thread count argument");
if(numthreads<1||numthreads> MAX_THREADS)
error("! bad value for thread count");
argc--;
argv++;
break;
/*5:*/
#line 87 "./twophase.w"

case'M':
if(sscanf(argv[1],"%lld",&phase2limit)!=1)
error("! bad argument to -M");
argc--;
argv++;
break;
case's':
if(sscanf(argv[1],"%d",&target_length)!=1)
error("! bad argument to -s");
if(target_length>=MAX_MOVES)
target_length= MAX_MOVES-1;
argc--;
argv++;
break;

/*:5*//*8:*/
#line 122 "./twophase.w"

case'W':skipwrite++;break;

/*:8*//*20:*/
#line 254 "./twophase.w"

case'a':
axesmask= atol(argv[1]);
argv++;
argc--;
break;

/*:20*/
#line 65 "./twophase.w"

default:
error("! bad argument");
}
}

/*:3*/
#line 28 "./twophase.w"

/*6:*/
#line 107 "./twophase.w"

if(phase2limit>=0xffffffffffffffLL&&target_length==0&&
verbose<=1)
error("! must specify -M, -s, or -v");

/*:6*//*9:*/
#line 130 "./twophase.w"

if(verbose)
cout<<BANNER<<endl<<flush;
phase1prune::init(skipwrite);
phase2prune::init(skipwrite);

/*:9*//*16:*/
#line 207 "./twophase.w"

pthread_mutex_init(&mutex,NULL);

/*:16*/
#line 29 "./twophase.w"

/*34:*/
#line 552 "./twophase.w"

pthread_t p_thread[MAX_THREADS];
for(int ti= 1;ti<numthreads;ti++)
pthread_create(&(p_thread[ti]),NULL,twophasesolver::worker,solvers+ti);
solvers[0].dowork();
for(int ti= 1;ti<numthreads;ti++)
pthread_join(p_thread[ti],0);/*:34*/
#line 30 "./twophase.w"

/*30:*/
#line 484 "./twophase.w"

if(missed_target)
cout<<"WARNING:  missed target on "<<missed_target
<<" sequences."<<endl;
phase1prune::check_integrity();
phase2prune::check_integrity();
cout<<"Solved "<<solved<<" sequences in "<<duration()
<<" seconds with "<<phase2total<<" probes."<<endl;
cout<<"Completed in "<<(walltime()-progstart)<<endl;
exit(0);

/*:30*/
#line 31 "./twophase.w"

}

/*:1*/
