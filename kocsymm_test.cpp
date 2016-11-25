/*62:*/
#line 1141 "./kocsymm.w"

#include "kocsymm.h"
#include <iostream> 
using namespace std;
int main(int argc,char*argv[]){
if(lrand48()==0)
srand48(getpid()+time(0));
kocsymm kc,kc2;
permcube pc,pc2;
cubepos cp,cp2;
/*56:*/
#line 1013 "./kocsymm.w"

{
cubepos cpi;
permcube pci(cpi);
kocsymm kci(cpi);
permcube pct;
kocsymm kct;
if(pct!=pci||kct!=kci)
error("! problem with default constructors");
if(permcube::c12_8[pc.et]!=0)
error("! bad mapping in 12->8");
}

/*:56*/
#line 1151 "./kocsymm.w"

/*57:*/
#line 1032 "./kocsymm.w"

for(int i= 0;i<100000;i++){
cp.randomize();
kocsymm kc(cp);
permcube pc(cp);
if(kc.epsymm!=pc.em)
error("! mismatch in edge middle occupancy");
kc.set_coset(cp2);
pc.set_perm(cp2);
if(cp!=cp2)
error("! mismatch in conversion and back");
}

/*:57*/
#line 1152 "./kocsymm.w"

/*58:*/
#line 1048 "./kocsymm.w"

for(int i= 0;i<1000;i++){
cp.randomize();
kocsymm kc(cp);
permcube pc(cp);
int mv= random_move_ext();
cp.movepc(mv);
cp2= cp;
kc.move(mv);
pc.move(mv);
kc.set_coset(cp2);
pc.set_perm(cp2);
if(cp!=cp2)
error("! mismatch in move test");
}

/*:58*/
#line 1153 "./kocsymm.w"

/*59:*/
#line 1066 "./kocsymm.w"

for(int i= 0;i<1000;i++){
cp.randomize();
kocsymm kc(cp);
for(int m= 1;m<KOCSYMM;m++){
cp.remap_into(m,cp2);
kocsymm kc2(cp2);
if(kc2<kc)
kc= kc2;
}
kocsymm kc3(cp);
kc3.canon_into(kc2);
if(kc2!=kc)
error("! canonicalization failuree");
}

/*:59*/
#line 1154 "./kocsymm.w"

/*60:*/
#line 1085 "./kocsymm.w"

int s= 0;
for(int c= 0;c<CORNERSYMM;c++){
int bits= kocsymm::cornersymm[c].minbits;
if(bits==1){
s+= EDGEOSYMM*EDGEPERM;
}else if(bits&1){
for(int eo= 0;eo<EDGEOSYMM;eo++)
for(int ep= 0;ep<EDGEPERM;ep++){
kocsymm kc(c,eo,ep);
kc.canon_into(kc2);
if(kc==kc2)
s++;
}
}
}
cout<<"Final sum is "<<s<<endl;
if(s!=138639780)
error("! bad total coset calculation");

/*:60*/
#line 1155 "./kocsymm.w"

/*61:*/
#line 1114 "./kocsymm.w"

int mvs[10000];
for(int i= 0;i<10000;i++)
mvs[i]= random_move();
duration();
for(int i= 0;i<10000;i++)
for(int j= 0;j<10000;j++)
kc.move(mvs[j]);
cout<<"Moving 100M kc took "<<duration()<<endl;
for(int i= 0;i<10000;i++)
for(int j= 0;j<10000;j++)
pc.move(mvs[j]);
cout<<"Moving 100M pc took "<<duration()<<endl;
for(int i= 0;i<10000;i++)
for(int j= 0;j<10000;j++)
cp.move(mvs[j]);
cout<<"Moving 100M cp (move) took "<<duration()<<endl;
for(int i= 0;i<10000;i++)
for(int j= 0;j<10000;j++)
cp.movepc(mvs[j]);
cout<<"Moving 100M cp (movepc) took "<<duration()<<endl;
if(cp<cp2&&pc<pc2&&kc<kc2)
cout<<"(Ignore this message.)"<<endl;


/*:61*/
#line 1156 "./kocsymm.w"

}/*:62*/
