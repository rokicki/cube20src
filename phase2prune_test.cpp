/*26:*/
#line 609 "./phase2prune.w"

#include "phase2prune.h"
#include <iostream> 
using namespace std;
char buf[4096];
int main(int argc,char*argv[]){
if(lrand48()==0)
srand48(time(0));
phase2prune::init(0);
phase2prune::check_integrity();
cubepos cp;
for(int i= 0;i<100000;i++){
char*tmp;
int mv= random_move();
if(kocsymm::in_Kociemba_group(mv)){
cp.movepc(mv);
}
int lookd= phase2prune::lookup(cp);
cout<<"Distance "<<lookd<<endl;
moveseq s= phase2prune::solve(cp);
cubepos cpt= cp;
for(unsigned int j= 0;j<s.size();j++)
cpt.movepc(s[j]);
cubepos::append_moveseq(tmp= buf,s);
cout<<"Solution length "<<s.size()<<" "<<buf<<endl;
if(cpt!=identity_cube)
error("! bad solve");
if((unsigned int)lookd> s.size())
error("! solution too short");
}
}/*:26*/
