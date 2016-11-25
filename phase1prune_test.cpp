/*23:*/
#line 475 "./phase1prune.w"

#include "phase1prune.h"
#include <iostream> 
using namespace std;
int main(int argc,char*argv[]){
if(lrand48()==0)
srand48(time(0));
phase1prune::init();
int t[10000];
for(int i= 0;i<10000;i++)
t[i]= random_move();
kocsymm kc;
/*24:*/
#line 495 "./phase1prune.w"

for(int i= 0;i<10000;i++){
kc.move(t[i]);
int d= phase1prune::lookup(kc);
for(int d2= d-1;d2<=d+2;d2++){
int nextmovemask= 0;
kocsymm kc2;
int dt= phase1prune::lookup(kc,d2,nextmovemask);
if(dt!=d)
error("mismatch on lookup");
int bc= 0;
for(int mv= 0;mv<NMOVES;mv++){
kc2= kc;
kc2.move(mv);
dt= phase1prune::lookup(kc2);
if(dt<d2)
bc|= (1<<mv);
}
if(nextmovemask!=bc)
error("! move mask error");
}
}

/*:24*/
#line 487 "./phase1prune.w"

/*25:*/
#line 520 "./phase1prune.w"

int sum= 0;
duration();
for(int i= 0;i<10000;i++)
for(int j= 0;j<10000;j++){
kc.move(t[j]);
sum+= phase1prune::lookup(kc);
}
cout<<"Did 100M basic lookups in "<<duration()<<" sum "<<sum<<endl;

/*:25*/
#line 488 "./phase1prune.w"

/*26:*/
#line 532 "./phase1prune.w"

int prev= 10;
int nextmovemask= 0;
for(int i= 0;i<10000;i++)
for(int j= 0;j<10000;j++){
kc.move(t[j]);
int r= phase1prune::lookup(kc,prev,nextmovemask);
sum+= r+nextmovemask;
prev= r;
}
cout<<"Did 100M extended lookups in "<<duration()<<" sum "<<sum<<endl;

/*:26*/
#line 489 "./phase1prune.w"

/*27:*/
#line 546 "./phase1prune.w"

for(int i= 0;i<1000;i++)
for(int j= 0;j<10000;j++){
kc.move(t[j]);
phase1prune::solve(kc);
}
cout<<"Did 10M solves in "<<duration()<<" sum "<<sum<<endl;/*:27*/
#line 490 "./phase1prune.w"

}

/*:23*/
