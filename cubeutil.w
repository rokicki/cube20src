@* Do miscellaneous operations on cube positions.

@(cubeutil.cpp@>=
#include "phase1prune.h"
#include <map>
#include <iostream>
#include <vector>
#include <map>
#include <cstdio>
using namespace std ;
char op = 0 ;
int n = 2000000000 ;
int uniq = 0 ;
map<cubepos, char *> world ;
int getparity(const cubepos &cp) {
   int r = 0 ;
   for (int i=0; i<8; i++)
      for (int j=i+1; j<8; j++)
         if (cubepos::corner_perm(cp.c[i]) > cubepos::corner_perm(cp.c[j]))
            r++ ;
   return r & 1 ;
}
map<long long, const char *> symmlookup ;
const char *base[] = {
"U1U3 //O_h",
"U2L2F2D2U2F2R2U2 //T_h",
"B1F1L1R1B3F3D3U3L1R1D1U1 //T",
"U1D3 //C_4h",
"U2 //C_4v",
"U1 //C_4",
"U1L3R3B2U3R2B1L2D3F2L3R3U3 //C_3v",
"B2D3U3R2B2U2F2D1U3 //C_3",
"U2R2D2U2R2D3U3 //C_2h_(a)",
"U2F2R2F2R2F2R2U2 //C_2h_(b)",
"B2F2U1R2B2F2R2U3L2R2 //C_2v_(a1)",
"B2F2U2 //C_2v_(a2)",
"D3B2R2B2D3U3F2L2F2U3 //C_2v_(b)",
"D2L2R2U3L2R2D3 //C_2_(a)",
"F2U1F2U1B2F2D3B2U3F2 //C_2_(b)",
"U2F2L2D2L2U2R2U2F2U2 //C_s_(a)",
"U2F2L2F2R2F2R2 //C_s_(b)",
"R2U2F2U2R2 //C_i",
"U1R1 //C_1",
"D2U2 //D_4h",
"U1D1 //D_4",
"U1L1D1U1L3D3U3R1B2U2B2L3R3U3 //D_3d",
"D1B1D1U2B2F2L2R2U3F1U1 //D_3",
"U3L2R2U2B2F2U3 //D_2h_(edge)",
"R2B2F2R2 //D_2h_(face)",
"B1F1L2R2B3F3U2 //D_2d_(edge)",
"U3L1R1B2F2L3R3U3 //D_2d_(face)",
"U3R2D2U2R2D2U3 //D_2_(edge)",
"R2B2L2B2D1U1B2L2F2R2D1U1 //D_2_(face)",
"B2R2D3U3B2R2D3U3B2R2 //S_4",
"B3D3U1L3R1B3F1U1 //S_6"
} ;
long long getsymm(const cubepos &cp) {
   cubepos cp2 ;
   long long r = 0 ;
   for (int m=0; m<48; m++) {
      cp.remap_into(m, cp2) ;
      if (cp2 == cp)
         r |= 1LL << m ;
   }
   return r ;
}
void parseposition(cubepos &cp, const char *s) {
  cp = identity_cube ;
  if (cp.parse_Singmaster(s) == 0)
     return ;
  cp = identity_cube ;
  const char *faces = "UuFfRrDdBbLlIiJjKk" ;
#ifdef AXIAL
  faces = "UuFfRrDdBbLlAaCcEeGgHhIiJjKkMm" ;
#endif
  while (*s) {
    while (*s && *s <= ' ')
       s++ ;
    if (*s == 0 || *s == '/')
       break ;
    const char *fp = strchr(faces, *s) ;
    if (fp == 0) {
       cerr << "Illegal face in " << *s << endl ;
       exit(10) ;
    }
    int f = (fp - faces) / 2 ;
    if (f * TWISTS > NMOVES)
       error("! bad face for this metric") ;
    s++ ;
    switch (*s) {
case '-':
case '3':
case '\'':
       cp.movepc(f*TWISTS+TWISTS-1) ;
       break ;
case '2':
#ifdef QUARTER
       cp.movepc(f*TWISTS) ;
       cp.movepc(f*TWISTS) ;
#else
       cp.movepc(f*TWISTS+1) ;
#endif
       break ;
case '+':
case '1':
       cp.movepc(f*TWISTS) ;
       break ;
default:
       cerr << "Illegal modifier in " << *s << endl ;
       exit(10) ;
    }
    s++ ;
  }
}
void initsymm() {
   for (int i=0; i<sizeof(base)/sizeof(base[0]); i++) {
      cubepos cp, cp2 ;
      parseposition(cp, base[i]) ;
      const char *p = index(base[i], '/') + 2 ;
      for (int m=0; m<48; m++) {
         cp.remap_into(m, cp2) ;
         long long r = getsymm(cp2) ;
         symmlookup[r] = p ;
      }
   }
}
char buf[1001] ;
const char *dict = 0 ;
int getpos(cubepos &cp, FILE *f=stdin) {
   if (fgets(buf, 1000, f) == 0)
      return 0 ;
   parseposition(cp, buf) ;
   return 1 ;
}
void obeyuniq(cubepos &cp2, const cubepos &cp) {
   if (uniq == 1) {
      cp2 = cp ;
   } else if (uniq == 48) {
      cp.canon_into48(cp2) ;
   } else if (uniq == 96) {
      cp.canon_into96(cp2) ;
   } else {
      error("! bad uniq value") ;
   }
}
int main(int argc, char *argv[]) {
   cubepos cp ;
   kocsymm kc ;
   while (argc > 1 && argv[1][0] == '-') {
      argc-- ;
      argv++ ;
      switch (argv[0][1]) {
case 'r':  // generate random positions
         op = 'r' ;
         n = atol(argv[1]) ;
         argc-- ;
         argv++ ;
         break ;
case 'p':  // print parity of input positions
         op = 'p' ;
         break ;
case 'k':  // generate random Kociemba coset positions
         op = 'k' ;
         n = atol(argv[1]) ;
         argc-- ;
         argv++ ;
         break ;
case 'u':  // uniquify; can be immediately followed by 1, 48, or 96
case 'd':  // difference:  only print positions from stdin that aren't in filea
case 'l':  // lookup:  find solutions based on database of positions
         uniq = 96 ;
         op = argv[0][1] ;
         if (argv[0][2])
            uniq = atol(argv[0]+2) ;
         if (op == 'd' || op == 'l') {
            argc-- ;
            argv++ ;
            dict = argv[0] ;
         }
         break ;
case 's': // symmetry
         op = 's' ;
         break ;
case 'S': // show Singmaster
         op = 'S' ;
         break ;
case 'x': // count correct perms
         op = 'x' ;
         break ;
case 'X': // expand to all orientations and inversions
         op = 'X' ;
         break ;
default:
         error("! bad argument") ;
      }
   }
   switch (op) {
case 'r':
      while (n--) {
         cp.randomize() ;
         cout << cp.Singmaster_string() << endl ;
      }
      break ;
case 'p':
      while (getpos(cp)) {
         int p = getparity(cp) ;
         cout << p << endl ;
      }
      break ;
case 'k':
      phase1prune::init() ;
      while (n--) {
         cp.randomize() ;
         kocsymm kc(cp) ;
         moveseq ms = phase1prune::solve(kc) ;
         moveseq msi = cubepos::invert_sequence(ms) ;
         cout << cubepos::moveseq_string(msi) << endl ;
      }
      break ;
case 's':
      initsymm() ;
      while (getpos(cp)) {
         long long r = getsymm(cp) ;
         const char *p = symmlookup[r] ;
         while (buf[0] && buf[strlen(buf)-1] < ' ')
            buf[strlen(buf)-1] = 0 ;
         cout << buf << " " << p ;
         cubepos cpi ;
         cp.invert_into(cpi) ;
         const char *ai = "-" ;
         if (cp == cpi)
            ai = "I" ;
         else {
            cubepos cp2 ;
            for (int m=1; m<48; m++) {
               cpi.remap_into(m, cp2) ;
               if (cp2 == cp) {
                  ai = "A" ;
                  break ;
               }
            }
         }
         cout << " " << ai ;
         int s = 0 ;
         while (r) {
            s++ ;
            r &= r-1 ;
         }
         cout << " " << s << endl ;
      }
      break ;
case 'S':
      while (getpos(cp)) {
         cout << cp.Singmaster_string() << endl ;
      }
      break ;
case 'x':
      while (getpos(cp)) {
         int badedge = 0 ;
         int badcorn = 0 ;
         for (int i=0; i<12; i++)
            if (cubepos::edge_ori(cp.e[i]) ||
                ((i & 4) != (4 & cubepos::edge_perm(cp.e[i]))))
               badedge++ ;
         for (int i=0; i<8; i++)
            if (cubepos::corner_ori(cp.c[i]))
               badcorn++ ;
         cout << (badedge+badcorn) << " " << badedge << " " << badcorn << endl ;
      }
      break ;
case 'X':
      while (getpos(cp)) {
         cubepos cpi ;
         for (int inv=0; inv<2; inv++) {
            for (int m=0; m<M; m++) {
               cp.remap_into(m, cpi) ;
               cout << cpi.Singmaster_string() << endl ;
            }
            cpi = cp ;
            cpi.invert_into(cp) ;
         }
      }
      break ;
case 'd':
case 'l':
      {
         FILE *f = fopen(dict, "r") ;
         if (f == 0)
            error("! can't open dict file") ;
         while (getpos(cp, f)) {
            cubepos cp2 ;
            obeyuniq(cp2, cp) ;
            if (world.find(cp2) == world.end())
               world[cp2] = strdup(buf) ;
         }
      }
      // now fall through down to u
case 'u':
      while (getpos(cp)) {
         cubepos cp2 ;
         obeyuniq(cp2, cp) ;
         map<cubepos, char *>::iterator it = world.find(cp2) ;
         if (it == world.end()) {
            if (op == 'u' || op == 'd') {
               char *p = strdup(buf) ;
               cout << p ;
               world[cp2] = p ;
            } else {
               cout << "No solution for " << buf ;
            }
         } else {
            if (op == 'l') {
               const char *pp = it->second ;
               moveseq ms = cubepos::parse_moveseq(pp) ;
               while (*pp && *pp <= ' ')
                  pp++ ;
               if (*pp)
                  error("! dictionary input not parsed as moveseq") ;
               moveseq ms2 = ms ;
               moveseq res ;
               int found = 0 ;
               for (int inv=0; found == 0 && inv<2; inv++) {
                  for (int m=0; found == 0 && m<M; m++) {
                     cp2 = identity_cube ;
                     for (int i=0; i<(int)ms2.size(); i++) {
                        ms2[i] = cubepos::move_map[m][ms[i]] ;
                        cp2.move(ms2[i]) ;
                     }
                     if (cp2 == cp) {
                        res = ms2 ;
                        found = 1 ;
                     }
                  }
                  ms = cubepos::invert_sequence(ms) ;
               }
               if (found == 0)
                  error("! solution not found") ;
               cout << cubepos::moveseq_string(res) << endl ;
            }
         }
      }
      break ;
default:
      error("! bad op") ;
   }
}
