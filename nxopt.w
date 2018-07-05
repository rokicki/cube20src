@* Introduction.
This program is a new optimal solver for Rubik's cube.  Much like
Korf's original solver, it uses iterated depth-first search with large
pruning tables.  It further develops these ideas with a number of
additional small tricks.  It currently supports twelve different sizes
of pruning tables that allow you trade off memory consumption and
disk space for solving speed.

@(nxopt.cpp@>=
#include "cubepos.h"
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <sys/time.h>
#include <map>
#include <set>
#include <vector>
#undef CHECK
#undef TRACE
#undef TESTTABLE
using namespace std ;
typedef unsigned long long ull ;
typedef long long ll ;
const int CORNERSYMM = 2187 ;
const int C12_4 = 495 ;
const int C8_4 = 70 ;
unsigned short entropy[3] = { 1, 2, 3 } ;
int base ;
int bidir = 1 ;
int allsols = 0 ;
int numthreads = 1 ;
int mindepth = 0 ;
int maxdepth = 40 ;
#ifndef QUARTER
int leftover = 0 ;
#endif
int symmetry = 0 ;
double order_mult[CANONSEQSTATES] ;
unsigned char c8_4crush[256] ;
unsigned char c8_4expand[C8_4] ;
const int CORNERUNIQ = 9930 ;
struct mind {
   unsigned char m, ind ;
} mindbuf[C8_4*139] ;
struct firstlev {
   unsigned int *p ;
   mind *v ;
   int levsiz, levoff ;
} firstlevs[CORNERSYMM] ;
int bc(int i) {
   int r = 0 ;
   while (i) {
      r++ ;
      i &= i-1 ;
   }
   return r ;
}
/**
 *   Faster corner coordinate computation.
 *
 *   This includes the |C8_4| multiplication, so requires 18 bits for the
 *   orientation.  Then the 4-of-8 is done in the next 8 bits, for a total
 *   of 26 bits.  So:
 *
 *   |fastcornercoords[i][j] = C8_4 * (j % 3) * 3^i (if i < 7) +
 *                            1 << (18 + i) (if j < 12)|
 *
 *   Note that we only look at the first *7*.
 */
int fastcornercoords[M][8][CUBIES] ;
int pow3[8] ;
void initfastcc() {
   pow3[0] = 1 ;
   for (int i=1; i<8; i++)
      pow3[i] = pow3[i-1] * 3 ;
   for (int m=0; m<M; m++) {
      int mprime = cubepos::invm[m] ;
      for (int ii=0; ii<7; ii++)
         for (int c=0; c<CUBIES; c++) {
            int c1 = cubepos::rot_corner[mprime][cubepos::corner_val(ii, 0)] ;
            int i = cubepos::corner_perm(c1) ;
            int c2 = cubepos::corner_ori_add(c, c1) ;
            int cc = cubepos::rot_corner[m][c2] ;
            fastcornercoords[m][i][c] = cubepos::corner_ori(cc) * pow3[ii] ;
            if (cubepos::corner_perm(cc) < 4)
               fastcornercoords[m][i][c] += 1<<(18+ii) ;
         }
   }
}
void getcornercoords(const cubepos &cp, int &cperm, int &cori, int m=0) {
   int s = 0 ;
   for (int i=0; i<8; i++)
      s += fastcornercoords[m][i][cp.c[i]] ;
   cperm = c8_4crush[s >> 18] ;
   cori = s & 0x3ffff ;
}
int getcornercoord(const cubepos &cp, int m=0) {
   int s = 0 ;
   for (int i=0; i<8; i++)
      s += fastcornercoords[m][i][cp.c[i]] ;
   return (s & 0x3ffff) * C8_4 + c8_4crush[s >> 18] ;
}
void setcornercoords(cubepos &cp, int cperm, int cori) {
   int orisum = 15 ;
   int lo = 0 ;
   int hi = 4 ;
   cperm = c8_4expand[cperm] ;
   for (int i=0; i<8; i++) {
      int mori = cori % 3 ;
      orisum -= mori ;
      if (i == 7)
         mori = orisum % 3 ;
      else
         cori /= 3 ;
      if ((cperm >> i) & 1) {
         cp.c[i] = cubepos::corner_val(lo, mori) ;
         lo++ ;
      } else {
         cp.c[i] = cubepos::corner_val(hi, mori) ;
         hi++ ;
      }
   }
}

@ In general, when looking up values, we use the following
array to accumulate information.

@(nxopt.cpp@>=
ull fastedges[M][12][CUBIES] ;

@ Now we turn our attention to the various ways we can do
edge orientation (E2).  There are only three ways:  EO1 uses
4 bits, EO2 uses 8 bits, and EO3 uses 11 bits.  We get the
low 11 bits of fastedges to work with to look up these values.
In the final index, the bits are located starting at least
significant bit nine (to save room for the low 9 bits giving
the edge occupancy).

@(nxopt.cpp@>=
const int POW3_6 = 729 ;
const int POW2_12 = 4096 ;
const int POW2_9 = 512 ;
#ifdef EO1
const int E2 = 16 ;
const int E2BITS = 4 ;
#else
#ifdef EO2
const int E2 = 256 ;
const int E2BITS = 8 ;
#else
#ifdef EO3
const int E2 = 2048 ;
const int E2BITS = 11 ;
#else
error "Please define one of EO1 EO2 or EO3" ;
#endif
#endif
#endif
#ifndef EO3
unsigned char eolo[POW3_6], eohi[POW3_6] ;
void initeorecur(int togo, int ind, int bits, int nbits) {
   if (togo == 0) {
      eolo[ind] = bits ;
      eohi[ind] = bits << (E2BITS - nbits) ;
   } else {
      initeorecur(togo-1, 3*ind, bits, nbits) ;
      initeorecur(togo-1, 3*ind+1, 2*bits, nbits+1) ;
      initeorecur(togo-1, 3*ind+2, 2*bits+1, nbits+1) ;
   }
}
#endif
void initeo() {
   for (int i=0; i<12; i++) {
      int ii = (i + 8) % 12 ;
      for (int c=0; c<CUBIES; c++) {
#ifdef EO1
         if (c & 8) { // middle cubie
            if (ii < 6)
               fastedges[0][i][c] = ((cubepos::edge_ori(c) + 1) * pow3[ii]) ;
            else
               fastedges[0][i][c] =
                            ((cubepos::edge_ori(c) + 1) * pow3[ii-6]) << 10 ;
         }
#endif
#ifdef EO2
         if (0 == (c & 8)) { // not middle cubie
            if (ii < 6)
               fastedges[0][i][c] = ((cubepos::edge_ori(c) + 1) * pow3[ii]) ;
            else
               fastedges[0][i][c] =
                            ((cubepos::edge_ori(c) + 1) * pow3[ii-6]) << 10 ;
         }
#endif
#ifdef EO3
         if (ii < 11)
            fastedges[0][i][c] = cubepos::edge_ori(c) << ii ;
#endif
      }
   }
#ifndef EO3
   initeorecur(6, 0, 0, 0) ;
#endif
}

@ Now we turn our attention to the various ways we can do
edge permutation (E1).  There are only four ways:  EP1 uses
only the middle edge occupancy; EP2 uses that and the
permutation of the middle edges; EP3 uses the middle and
bot edge occupancy, and EP4 uses the middle and bot edge
occupancy, as well as the permutation of the middle edges.
We get the 34 bits starting at least significant bit 10
to work with. For middle edge occupany, we use two six-bit
fields.  For middle and bot edge occupancy, we use two
10-bit fields.  For middle edge permutation, we use
one 14-bit field.  We always place the occupany bits
starting at least significant bit 20, and the permutation
bits starting at least significant bit 40.
The E1 value takes into account the 495->512 expansion.  

@(nxopt.cpp@>=
#ifdef EP1
const int E1 = POW2_9 ;
#else
#ifdef EP2
const int E1 = POW2_9 * 24 ;
#define USEEPERM
#else
#ifdef EP3
const int E1 = POW2_9 * C8_4 ;
#define USEEBOT
#else
#ifdef EP4
const int E1 = POW2_9 * C8_4 * 24 ;
#define USEEPERM
#define USEEBOT
#else
error "Please define one of EP1, EP2, EP3, or EP4."
#endif
#endif
#endif
#endif
char fmepermi[C12_4] ; // just during init
#ifdef USEEPERM
unsigned char fmeperm[12*12*12*6] ;
#endif
#ifdef USEEBOT
unsigned short emphi[POW3_6] ;
unsigned short emplo[POW3_6] ;
unsigned char ebphi[POW3_6] ;
unsigned char ebplo[POW3_6] ;
#else
unsigned short ephi[64] ;
unsigned char eplo[64] ;
#endif
unsigned short bitstoindex[POW2_12] ;
unsigned short c12_4expand[C12_4] ;
unsigned short pack11_9[2048] ;
unsigned char midperm[24][4] ;
void initfmeperm(int bits, int ind, int togo) {
   if (togo == 0) {
      int perm = --fmepermi[bitstoindex[bits]] ;
#ifdef USEEPERM
      fmeperm[ind>>1] = perm ;
#endif
      if (bits == 15) {
         midperm[perm][ind/1728] = 14 ;
         midperm[perm][ind/144%12] = 12 ;
         midperm[perm][ind/12%12] = 10 ;
         midperm[perm][ind%12] = 8 ;
      }
   } else {
      for (int i=0; i<12; i++)
         if (((bits >> i) & 1) == 0)
            initfmeperm(bits|(1<<i), 12*ind+i, togo-1) ;
   }
}
void initeprecur(int togo, int ind, int bits, int nbits, int bits2, int nbits2) {
   if (togo == 0) {
      if (nbits <= 4 && nbits2 <= 4) {
#ifdef USEEBOT
         emplo[ind] = bitstoindex[bits] ;
         emphi[ind] = bitstoindex[(bits<<6)+(1<<(4-nbits))-1] ;
         ebplo[ind] = bitstoindex[bits2] ;
         ebphi[ind] = bitstoindex[(bits2<<(nbits+2))+(1<<(4-nbits2))-1] ;
#else
         ephi[bits] = bitstoindex[(bits<<6)+(1<<(4-nbits))-1] ;
         eplo[bits] = bitstoindex[bits] ;
#endif
      }
   } else {
      initeprecur(togo-1, 3*ind, bits*2+1, nbits+1, bits2, nbits2) ;
      initeprecur(togo-1, 3*ind+1, bits*2, nbits, bits2*2+1, nbits2+1) ;
      initeprecur(togo-1, 3*ind+2, bits*2, nbits, bits2*2, nbits2) ;
   }
}
void initep() {
   for (int i=0; i<12; i++) {
      int ii = (i + 8) % 12 ;
      for (int c=0; c<CUBIES; c++) {
#ifdef USEEBOT
         int v = (2 + (c >> 3)) % 3 ;
         if (ii < 6)
            fastedges[0][i][c] += ((ull)(v * pow3[ii])) << 20 ;
         else
            fastedges[0][i][c] += ((ull)(v * pow3[ii-6])) << 30 ;
#else
         fastedges[0][i][c] += ((ull)((c >> 3) & 1)) << (20 + ii) ;
#endif
#ifdef USEEPERM
         if (c & 8) {
            int ceperm = cubepos::edge_perm(c)-4 ;
            fastedges[0][i][c] +=
                  ((ull)(((ii * pow3[ceperm]) << (2 * ceperm)) >> 1)) << 40 ;
         }
#endif
   if (fastedges[0][i][c] > 1LL << 54)
      error("! bad fastedges") ;
      }
   }
   int bitval[13] ;
   for (int i=0; i<13; i++)
      bitval[i] = 0 ;
   for (int i=0; i<POW2_12; i++) {
      int nbits = bc(i) ;
      if (nbits == 4)
         c12_4expand[bitval[nbits]] = i ;
      bitstoindex[i] = bitval[nbits]++ ;
   }
   int at = 0 ;
   for (int i=0; i<POW2_12; i++)
      if (bc(i) == 4)
         pack11_9[i&2047] = at++ ;
   initeprecur(6, 0, 0, 0, 0, 0) ;
   for (int i=0; i<C12_4; i++)
      fmepermi[i] = 24 ;
   initfmeperm(0, 0, 4) ;
}

@ Now we have the setters and getters.

@(nxopt.cpp@>=
#ifdef HALF
static const char *metric = "HALF" ;
#ifdef EP1
#ifdef EO1
#define BASE 7
#define DATFILE "nxopth11b.dat"
#endif
#ifdef EO2
#define BASE 8
#define DATFILE "nxopth21b.dat"
#endif
#ifdef EO3
#define BASE 9
#define DATFILE "nxopth31b.dat"
#endif
#endif
#ifdef EP2
#ifdef EO1
#define BASE 8
#define DATFILE "nxopth12b.dat"
#endif
#ifdef EO2
#define BASE 9
#define DATFILE "nxopth22b.dat"
#endif
#ifdef EO3
#define BASE 10
#define DATFILE "nxopth32b.dat"
#endif
#endif
#ifdef EP3
#ifdef EO1
#define BASE 8
#define DATFILE "nxopth13b.dat"
#endif
#ifdef EO2
#define BASE 9
#define DATFILE "nxopth23b.dat"
#endif
#ifdef EO3
#define BASE 10
#define DATFILE "nxopth33b.dat"
#endif
#endif
#ifdef EP4
#ifdef EO1
#define BASE 9
#define DATFILE "nxopth14b.dat"
#endif
#ifdef EO2
#define BASE 11 // ??
#define DATFILE "nxopth24b.dat"
#endif
#ifdef EO3
#define BASE 11 // ??
#define DATFILE "nxopth34b.dat"
#endif
#endif
#endif
#ifdef QUARTER
static const char *metric = "QUARTER" ;
#ifdef EP1
#ifdef EO1
#define BASE 8
#define DATFILE "nxoptq11b.dat"
#endif
#ifdef EO2
#define BASE 9
#define DATFILE "nxoptq21b.dat"
#endif
#ifdef EO3
#define BASE 10
#define DATFILE "nxoptq31b.dat"
#endif
#endif
#ifdef EP2
#ifdef EO1
#define BASE 9
#define DATFILE "nxoptq12b.dat"
#endif
#ifdef EO2
#define BASE 10
#define DATFILE "nxoptq22b.dat"
#endif
#ifdef EO3
#define BASE 11
#define DATFILE "nxoptq32b.dat"
#endif
#endif
#ifdef EP3
#ifdef EO1
#define BASE 10
#define DATFILE "nxoptq13b.dat"
#endif
#ifdef EO2
#define BASE 11
#define DATFILE "nxoptq23b.dat"
#endif
#ifdef EO3
#define BASE 12
#define DATFILE "nxoptq33b.dat"
#endif
#endif
#ifdef EP4
#ifdef EO1
#define BASE 11
#define DATFILE "nxoptq14b.dat"
#endif
#ifdef EO2
#define BASE 13 // ??
#define DATFILE "nxoptq24b.dat"
#endif
#ifdef EO3
#define BASE 13 // ??
#define DATFILE "nxoptq34b.dat"
#endif
#endif
#endif
#ifdef AXIAL
static const char *metric = "AXIAL" ;
#ifdef EP1
#ifdef EO1
#define BASE 5
#define DATFILE "nxopta11b.dat"
#endif
#ifdef EO2
#define BASE 6
#define DATFILE "nxopta21b.dat"
#endif
#ifdef EO3
#define BASE 6
#define DATFILE "nxopta31b.dat"
#endif
#endif
#ifdef EP2
#ifdef EO1
#define BASE 6
#define DATFILE "nxopta12b.dat"
#endif
#ifdef EO2
#define BASE 6
#define DATFILE "nxopta22b.dat"
#endif
#ifdef EO3
#define BASE 7
#define DATFILE "nxopta32b.dat"
#endif
#endif
#ifdef EP3
#ifdef EO1
#define BASE 6
#define DATFILE "nxopta13b.dat"
#endif
#ifdef EO2
#define BASE 7
#define DATFILE "nxopta23b.dat"
#endif
#ifdef EO3
#define BASE 7
#define DATFILE "nxopta33b.dat"
#endif
#endif
#ifdef EP4
#ifdef EO1
#define BASE 7
#define DATFILE "nxopta14b.dat"
#endif
#ifdef EO2
#define BASE 8 // ??
#define DATFILE "nxopta24b.dat"
#endif
#ifdef EO3
#define BASE 8 // ??
#define DATFILE "nxopta34b.dat"
#endif
#endif
#endif
#ifdef SLICE
static const char *metric = "SLICE" ;
#ifdef EP1
#ifdef EO1
#define BASE 6
#define DATFILE "nxopts11b.dat"
#endif
#ifdef EO2
#define BASE 7
#define DATFILE "nxopts21b.dat"
#endif
#ifdef EO3
#define BASE 8
#define DATFILE "nxopts31b.dat"
#endif
#endif
#ifdef EP2
#ifdef EO1
#define BASE 7
#define DATFILE "nxopts12b.dat"
#endif
#ifdef EO2
#define BASE 8
#define DATFILE "nxopts22b.dat"
#endif
#ifdef EO3
#define BASE 9
#define DATFILE "nxopts32b.dat"
#endif
#endif
#ifdef EP3
#ifdef EO1
#define BASE 7
#define DATFILE "nxopts13b.dat"
#endif
#ifdef EO2
#define BASE 8
#define DATFILE "nxopts23b.dat"
#endif
#ifdef EO3
#define BASE 9
#define DATFILE "nxopts33b.dat"
#endif
#endif
#ifdef EP4
#ifdef EO1
#define BASE 8
#define DATFILE "nxopts14b.dat"
#endif
#ifdef EO2
#define BASE 9 // ??
#define DATFILE "nxopts24b.dat"
#endif
#ifdef EO3
#define BASE 9 // ??
#define DATFILE "nxopts34b.dat"
#endif
#endif
#endif
#ifdef OVERRIDEBASE
#undef BASE
#define BASE OVERRIDEBASE
#endif
int getedgecoord(const cubepos &cp, int m=0) {
   ull s = 0 ;
   for (int i=0; i<12; i++)
      s += fastedges[m][i][cp.e[i]] ;
#ifdef EO3
   int eoc = s & 0x7ff ;
#else
   int eoc = eolo[s&0x3ff] + eohi[(s>>10)&0x3ff] ;
#endif
#ifdef USEEBOT
   int ep1 = (s>>20)&0x3ff ;
   int ep2 = (s>>30)&0x3ff ;
   int epmc = emplo[ep1] + emphi[ep2] ;
   int epbotc = ebplo[ep1] + ebphi[ep2] ;
#else
   int epmc = pack11_9[(s>>20)&0x7ff] ;
   int epbotc = 0 ;
#endif
#ifdef USEEPERM
   epbotc = 24 * epbotc + fmeperm[s>>40] ;
#endif
   int epadd = 2 * ((epmc + (epmc >> 5) + 65) >> 6) ;
   return (epbotc << (E2BITS + 9)) + (eoc << 9) + epmc + epadd ;
}
/*
 *   This is not particularly speed-sensitive.
 */
void setedgecoord(cubepos &cp, int ecoord) {
   int emocc = ecoord & 0x1ff ;
   emocc -= 2 + 2 * (emocc >> 6) ;
#ifndef EO3
   int paritybit = 0 ;
#endif
#ifdef EO1
   int eoperm = (ecoord >> 9) & 0xf ;
   if (bc(eoperm) & 1)
      paritybit = 1 ;
#endif
#ifdef EO2
   int eoperm = (ecoord >> 9) & 0xff ;
   if (bc(eoperm) & 1)
      paritybit = 1 ;
#endif
#ifdef EO3
   int eoperm = (ecoord >> 9) & 0x7ff ;
   if (bc(eoperm) & 1)
      eoperm |= 0x800 ;
#endif
   int ebocc = ecoord >> (E2BITS + 9) ;
#ifdef USEEPERM
   int emperm = ebocc % 24 ;
#ifdef USEEBOT
   ebocc /= 24 ;
#endif
#else
   int emperm = 0 ;
#endif
   unsigned char *mid = midperm[emperm] ;
   emocc = c12_4expand[emocc] ;
   ebocc = c8_4expand[ebocc] ;
 if (bc(emocc) != 4) error("! bad emocc") ;
 if (bc(ebocc) != 4) error("! bad ebocc") ;
   static unsigned char topedges[4] = { 0, 2, 4, 6 } ;
   static unsigned char botedges[4] = { 16, 18, 20, 22 } ;
   unsigned char *top = topedges ;
   unsigned char *bot = botedges ;
   for (int i=0; i<12; i++, emocc >>= 1) {
      int ii = (i + 4) % 12 ;
      if (emocc & 1) {
         cp.e[ii] = *mid++ ;
#ifndef EO2
         cp.e[ii] |= eoperm & 1 ;
         eoperm >>= 1 ;
#else
         cp.e[ii] |= paritybit ;
         paritybit = 0 ;
#endif
      } else {
         if (ebocc & 1) {
            cp.e[ii] = *bot++ ;
         } else {
            cp.e[ii] = *top++ ;
         }
         ebocc >>= 1 ;
#ifndef EO1
         cp.e[ii] |= eoperm & 1 ;
         eoperm >>= 1 ;
#else
         cp.e[ii] |= paritybit ;
         paritybit = 0 ;
#endif
      }
   }
}

@ Expand fastedges for other mappings.

@(nxopt.cpp@>=
void initfillout() {
   for (int m=1; m<M; m++) {
      int mprime = cubepos::invm[m] ;
      for (int i=0; i<12; i++)
         for (int c=0; c<CUBIES; c++) {
            int c1 = cubepos::rot_edge[mprime][i*2] ;
            int i2 = cubepos::edge_perm(c1) ;
            int c2 = cubepos::edge_ori_add(c, c1) ;
            fastedges[mprime][i][cubepos::rot_edge[m][c2]] = fastedges[0][i2][c] ;
         }
   }
}

@ Finally we have more common code.

@(nxopt.cpp@>=
const int EDGECOOR = E1 * E2 ;
const size_t BIGMEMSIZE = (EDGECOOR+3)/4 ;
const size_t MEMOFFSET = (BIGMEMSIZE+3)/4 ;
const size_t MEMSIZE = 4 * MEMOFFSET ;
void initmem() {
   mind *mindp = mindbuf ;
   int totsiz = 0 ;
   for (int j=0; j<CORNERSYMM; j++) {
      mind mindt[C8_4] ;
      memset(mindt, 0, sizeof(mindt)) ;
      int siz = 0 ;
      int mj = j ;
      for (int i=0; i<C8_4; i++) {
         int c = j*C8_4+i ;
         int rm = 0 ;
         for (int k=1; k<16; k++) {
            cubepos cp, cp2 ;
            int ii, jj ;
            setcornercoords(cp, i, j) ;
            cp.remap_into(k, cp2) ;
            getcornercoords(cp2, ii, jj) ;
            int tc = jj*C8_4+ii ;
            if (tc < c) {
               c = tc ;
               rm = k ;
               mj = jj ;
            }
         }
         mindt[i].m = rm ;
         if (c / C8_4 == j) {
            if (c % C8_4 == i) {
               mindt[i].ind = siz++ ;
            } else {
               mindt[i].ind = mindt[c % C8_4].ind ;
            }
         } else {
            mindt[i].ind = firstlevs[c/C8_4].v[c % C8_4].ind ;
            siz = -1 ;
         }
      }
      firstlevs[j].levsiz = siz ;
      for (int i=0; i<j; i++)
         if (memcmp(mindt, firstlevs[i].v, sizeof(mindt)) == 0) {
            firstlevs[j].v = firstlevs[i].v ;
            break ;
         }
      if (firstlevs[j].v == 0) {
         firstlevs[j].v = mindp ;
         mindp += C8_4 ;
         memcpy(firstlevs[j].v, mindt, sizeof(mindt)) ;
      }
      if (siz > 0) {
         firstlevs[j].p = (unsigned int *)malloc(MEMSIZE * (size_t)siz) ;
         firstlevs[j].levoff = totsiz ;
         totsiz += siz ;
      } else {
         firstlevs[j].p = firstlevs[mj].p ;
         firstlevs[j].levoff = firstlevs[mj].levoff ;
      }
   }
   if (totsiz != CORNERUNIQ) {
      error("! unexpected totsiz") ;
   }
}
const int SIZE = 1000000 ;
unsigned char tomove[SIZE] ;
unsigned char tom[SIZE] ;
void test() {
   cubepos cp, cp2, cp3 ;
   for (int trial=0; trial<1000; trial++) {
      int c = getedgecoord(cp) ;
      setedgecoord(cp2, c) ;
      int c2 = getedgecoord(cp2) ;
      if (c != c2)
         error("! bad match") ;
      int mv = random_move() ;
      cp.movepc(mv) ;
      while (1) {
         c = (int)(drand48()*E1*E2) ;
         if ((c & 63) >= 2 && (c & 511) < 511)
            break ;
      }
      setedgecoord(cp2, c) ;
      c2 = getedgecoord(cp2) ;
      if (c != c2) {
         error("! bad match 2") ;
      }
      int m = (int)(16*drand48()) ;
      c = getedgecoord(cp2, m) ;
      cp2.remap_into(m, cp3) ;
      c2 = getedgecoord(cp3, 0) ;
      if (c != c2)
         cout << "! botch in remap at " << m << endl ;
   }
   for (int i=0; i<SIZE; i++) {
      tomove[i] = random_move() ;
      tom[i] = (int)(16*drand48()) ;
   }
   duration() ;
   int s = 0 ;
   for (int j=0; j<10; j++)
   for (int i=0; i<SIZE; i++) {
      int cperm, cori ;
      getcornercoords(cp, cperm, cori) ;
      s += cperm + cori ;
      cp.move(tomove[i]) ;
   }
   cout << "Corner in " << duration() << " " << s << endl ;
   for (int j=0; j<10; j++)
   for (int i=0; i<SIZE; i++) {
      s += getedgecoord(cp, tom[i]) ;
      cp.move(tomove[i]) ;
   }
   cout << "Edge in " << duration() << " " << s << endl ;
}
unsigned char skipata[NMOVES] ;
long long expandm[8] ;
void init() {
   int j = 0 ;
   for (int i=0; i<256; i++) 
      if (bc(i) == 4) {
          c8_4crush[i] = j ;
          c8_4expand[j] = i ;
          c8_4crush[i & 127] = j ; // don't require last one
          j++ ;
      }
   for (int mv=0; mv<NMOVES; mv++)
      skipata[mv] = mv / TWISTS % 3 ;
   for (int i=0; i<8; i++)
      for (int mv=0; mv<NMOVES; mv++)
         if (((i >> skipata[mv]) & 1) == 0)
            expandm[i] |= 1LL << mv ;
   initfastcc() ;
   initeo() ;
   initep() ;
   initfillout() ;
// test() ;
   initmem() ;
}
int popcount64(long long v) {
   return __builtin_popcountll(v) ;
}
struct efast {
   int base ;
   int *bits ;
} emove[NMOVES][E1], emap[16][E1] ;
map<vector<int>, int*> e2offmap ;
int *finde2bits(int *bits) {
   vector<int> key ;
   for (int bi=0; bi<E2BITS; bi++)
      key.push_back(bits[bi]) ;
   if (e2offmap.find(key) == e2offmap.end()) {
      int *bits2 = (int *)calloc(E2, sizeof(int)) ;
      for (int bi=0; bi<E2BITS; bi++)
         for (int i=1<<bi; i<E2; i=(i+1)|(1<<bi))
            bits2[i] ^= bits[bi] ;
      e2offmap[key] = bits2 ;
   }
   return e2offmap[key] ;
}
void calcecoords() {
   cubepos cp, cp2, cp3, cp4 ;
   int bits[11] ;
   for (int ep=0; ep<E1; ep++) {
      if ((ep & 63) < 2 || (ep & 511) == 511)
         continue ;
      int baseep = (ep & 511) | ((ep & ~511) << E2BITS) ;
      for (int mv=0; mv<NMOVES; mv++) {
         setedgecoord(cp, baseep) ;
         cp.movepc(mv) ;
         int dec = getedgecoord(cp) ;
         emove[mv][ep].base = dec ;
         for (int bi=0, eo=1; eo<E2; eo += eo, bi++) {
            setedgecoord(cp3, baseep + (eo << 9)) ;
            cp3.movepc(mv) ;
            bits[bi] = dec ^ getedgecoord(cp3) ;
         }
         emove[mv][ep].bits = finde2bits(bits) ;
      }
      for (int m=0; m<16; m++) {
         setedgecoord(cp, baseep) ;
         cp.remap_into(m, cp2) ;
         int dec = getedgecoord(cp2) ;
         emap[m][ep].base = dec ;
         for (int bi=0, eo=1; eo<E2; eo += eo, bi++) {
            setedgecoord(cp3, baseep + (eo << 9)) ;
            cp3.remap_into(m, cp4) ;
            bits[bi] = dec ^ getedgecoord(cp4) ;
         }
         emap[m][ep].bits = finde2bits(bits) ;
      }
   }
}
long long have = 0, smhave = 0 ;
int globald ;
void dorow(unsigned int *srcp, long long &local_have, long long &local_smhave,
           unsigned int *dstp, int d3, int mv, int m) {
   int ds = (d3 + 1) % 3 ;
   int ec = 0 ;
   cubepos cp, cp2 ;
   efast *emovemv = emove[mv] ;
   efast *emapm = emap[m] ;
   for (int ep=0; ep<E1; ep += 512) {
      for (int eo=0; eo<E2; eo++) {
         for (int epm=0; epm<511; epm++, ec++) {
            if ((epm & 63) == 0 && (srcp[ec>>4] & 15) >= (unsigned int)globald) {
               epm += 63 ;
               ec += 63 ;
               continue ;
            }
            if ((epm & 15) == 0 && srcp[ec>>4] == 0xffffffff) {
               epm += 15 ;
               ec += 15 ;
               continue ;
            }
            if ((epm & 63) < 2)
               continue ;
            if ((int)((srcp[ec>>4] >> (2*(ec & 15))) & 3) == d3) {
               int e2 = emovemv[ep+epm].base ^ emovemv[ep+epm].bits[eo] ;
               int ep2 = (e2 & 511) + ((e2 >> E2BITS) & ~511) ;
               int eo2 = (e2 >> 9) & (E2 - 1) ;
               int dec = emapm[ep2].base ^ emapm[ep2].bits[eo2] ;
               if (((dstp[dec>>4] >> (2*(dec & 15))) & 3) == 3) {
                  local_have++ ;
                  dstp[dec>>4] -= (3 - ds) << (2*(dec & 15)) ;
                  if ((dstp[(dec&~63)>>4] & 15) == 15) {
                     dstp[(dec&~63)>>4] -= 15-globald ;
                     local_smhave++ ;
                  }
               }
            }
         }
         if (ec & 15)
            ec++ ; // take care of unused 511
      }
   }
   if (ec != E1 * E2)
      error("! oops 12") ;
}
void writetab() {
   duration() ;
   FILE *f = fopen(DATFILE, "wb") ;
   if (f == 0)
      error("! cannot write file") ;
   fputc('N', f) ;
   fputc('X', f) ;
   fputc(*metric, f) ;
   fputc(base, f) ;
   for (int i=0; i<CORNERSYMM; i++)
      for (int j=0; j<firstlevs[i].levsiz; j++)
         if ((int)fwrite(firstlevs[i].p+j*MEMOFFSET, 1, BIGMEMSIZE, f)
                                                           != BIGMEMSIZE)
            error("! error writing table") ;
   fclose(f) ;
   cout << "Table written in " << duration() << endl << flush ;
}
char colocks[CORNERUNIQ] ;
char codone[CORNERUNIQ] ;
int cocori[CORNERUNIQ] ;
int coperm[CORNERUNIQ] ;
int neighbors[CORNERUNIQ][NMOVES+1] ;
int coorder[CORNERUNIQ] ;
int genat ;
void calcodata() {
   cubepos cp, cp2 ;
   for (int cc=0; cc<CORNERSYMM; cc++) {
      firstlev &fl = firstlevs[cc] ;
      for (int j=0; j<C8_4; j++) {
         if (fl.v[j].m == 0) {
            int src = fl.v[j].ind + fl.levoff ;
            coorder[src] = src ;
            cocori[src] = cc ;
            coperm[src] = j ;
            neighbors[src][0] = src ;
            for (int mv=0; mv<NMOVES; mv++) {
               int dori, dperm ;
               setcornercoords(cp, j, cc) ;
               cp.movepc(mv) ;
               getcornercoords(cp, dperm, dori) ;
               firstlev &fld = firstlevs[dori] ;
               cp.remap_into(fld.v[dperm].m, cp2) ;
               neighbors[src][mv+1] = fld.v[dperm].ind + fld.levoff ;
            }
         }
      }
   }
   random_shuffle(coorder, coorder+CORNERUNIQ) ;
}
void startgenthreads() {
   memset(colocks, 0, sizeof(colocks)) ;
   memset(codone, 0, sizeof(codone)) ;
   genat = 0 ;
}
int getgenwork() {
   int r = -1 ;
   get_global_lock() ;
   for (int ii=genat; ii<CORNERUNIQ; ii++) {
      int i = coorder[ii] ;
      if (codone[i])
         continue ;
      int okay = 1 ;
      for (int j=0; j<NMOVES+1; j++)
         if (colocks[neighbors[i][j]])
            okay = 0 ;
      if (okay) {
         for (int j=0; j<NMOVES+1; j++)
            colocks[neighbors[i][j]] = 1 ;
         r = i ;
         break ;
      }
   }
   release_global_lock() ;
   return r ;
}
void finishgenwork(int i, long long local_have, long long local_smhave) {
   get_global_lock() ;
   have += local_have ;
   smhave += local_smhave ;
   while (genat < CORNERUNIQ && codone[coorder[genat]])
      genat++ ;
   for (int j=0; j<NMOVES+1; j++)
      colocks[neighbors[i][j]] = 0 ;
   codone[i] = 1 ;
   release_global_lock() ;
}
void dogenouter(int d3) {
   cubepos cp, cp2 ;
   while (1) {
      int r = getgenwork() ;
      if (r < 0)
         break ;
      long long local_have = 0, local_smhave = 0 ;
      int cc = cocori[r] ;
      int j = coperm[r] ;
      firstlev &fl = firstlevs[cc] ;
      unsigned int *srcp = fl.p + fl.v[j].ind*MEMOFFSET ;
      for (int mv=0; mv<NMOVES; mv++) {
         int dori, dperm ;
         setcornercoords(cp, j, cc) ;
         cp.movepc(mv) ;
         getcornercoords(cp, dperm, dori) ;
         cp.remap_into(firstlevs[dori].v[dperm].m, cp2) ;
         firstlev &fld = firstlevs[dori] ;
         int goalc = getcornercoord(cp2) ;
         for (int m=0; m<16; m++) {
            if (goalc == getcornercoord(cp, m))
               dorow(srcp, local_have, local_smhave,
                    fld.p+fld.v[dperm].ind*MEMOFFSET, d3, mv, m) ;
         }
      }
      finishgenwork(r, local_have, local_smhave) ;
   }
}
void *gworker(void *s) {
   int d3 = *(int *)s ;
   dogenouter(d3) ;
   return 0 ;
}
void generatetab() {
   calcecoords() ;
   calcodata() ;
   for (int i=0; i<CORNERSYMM; i++)
      for (int j=0; j<firstlevs[i].levsiz; j++)
         memset(firstlevs[i].p+j*MEMOFFSET, -1, BIGMEMSIZE) ;
   cubepos cp = identity_cube ;
   cubepos cp2 ;
   if (getcornercoord(cp) != 0)
      error("! bad initial state 2") ;
   if (getedgecoord(cp) != 2) {
      cout << "Initial state " << getedgecoord(cp) << endl ;
      error("! bad initial state 3") ;
   }
   firstlevs[0].p[0] &= ~0x3f ;
   ull didhave = 0 ;
   have = 1 ;
   smhave = 1 ;
   for (int d=0; d <= base+1; d++) {
      globald = d + 1 ;
      duration() ;
      didhave = have ;
      int d3 = d % 3 ;
      if (d >= base) {
         if (d == base) {
            for (int scc=0; scc<CORNERSYMM; scc++)
               for (int j=0; j<firstlevs[scc].levsiz; j++) {
                  unsigned int *srcp = firstlevs[scc].p+j*MEMOFFSET ;
                  for (int i=0; i<(int)BIGMEMSIZE; i += 16) {
                     unsigned int w = *srcp ;
                     unsigned int t = ((w & (w >> 1)) & 0x55555555) * 3 ;
                     *srcp++ = (t & ~0xf) + (w & 0xf) ;
                     w = *srcp ;
                     *srcp++ = ((w & (w >> 1)) & 0x55555555) * 3 ;
                     w = *srcp ;
                     *srcp++ = ((w & (w >> 1)) & 0x55555555) * 3 ;
                     w = *srcp ;
                     *srcp++ = ((w & (w >> 1)) & 0x55555555) * 3 ;
                  }
               }
         }
         d3 = d - base ;
      }
      startgenthreads() ;
#ifdef THREADS
      pthread_t p_thread[MAX_THREADS] ;
      for (int ti=0; ti<numthreads; ti++) {
         pthread_create(&(p_thread[ti]), NULL, gworker, (void*)&d3) ;
      }
      for (int ti=0; ti<numthreads; ti++) {
         pthread_join(p_thread[ti], 0) ;
      }
#else
      dogenouter(d3) ;
#endif
      cout << "At " << (d+1) << " have " << smhave << " " << have << " in "
           << duration() << endl ;
      if (have == (long long)didhave)
         break ;
   }
#ifndef TESTTABLE
   writetab() ;
#endif
}
int readtab() {
   duration() ;
   FILE *f = fopen(DATFILE, "rb") ;
   if (f == 0)
      return 0 ;
   if (fgetc(f) != 'N' || fgetc(f) != 'X')
      error("! bad datafile") ;
   if (fgetc(f) != *metric)
      error("! bad metric in datafile") ;
   int b = fgetc(f) ;
   if (b != base)
      printf("Overriding base from %d to %d\n", base, b) ;
   base = b ;
   cout << "Load " << DATFILE << " " << flush ;
   int goalperc = 2 ;
   int cnt = 0 ;
   for (int i=0; i<CORNERSYMM; i++) {
      for (int j=0; j<firstlevs[i].levsiz; j++) {
         if ((int)fread(firstlevs[i].p+j*MEMOFFSET, 1, BIGMEMSIZE, f)
                                                           != BIGMEMSIZE)
            error("! error reading table 1") ;
         cnt++ ;
         if (cnt * 100 >= goalperc * CORNERUNIQ) {
            if (goalperc % 10 == 0)
               cout << goalperc << flush ;
            else
               cout << "." << flush ;
            goalperc += 2 ;
         }
      }
   }
   fclose(f) ;
   cout << " read in " << duration() << endl << flush ;
   return 1 ;
}

@ More code.

@(nxopt.cpp@>=
void showmove(int mv) {
   char buf[20] ;
   char *p = buf ;
   cubepos::append_move(p, mv) ;
   cout << " " << buf ;
}
int getwork(cubepos &cp) ;
struct solution {
   int seq, length ;
   long long evals, probes ;
   double duration ;
   vector<char> moves ;
   void report() {
      if (length >= 0) {
         cout << "Solved " << seq << " in " << length ;
      } else {
         cout << "Did not solve " << seq << " in" ;
      }
      cout << " probes " << probes << " evals " << evals << " time " <<
          duration << endl << flush ;
      for (int i=0; i<moves.size(); i++) {
         showmove(moves[i]) ;
         if (length < 2 || (i + 1) % length == 0)
            cout << endl << flush ;
      }
   }
} ;
int outseq = 1 ;
vector<solution> pending_solutions ;
void report(solution &sol) {
   if (sol.seq == outseq) {
      sol.report() ;
      outseq++ ;
      while (1) {
         int ind = -1 ;
         for (int i=0; i<(int)pending_solutions.size(); i++)
            if (pending_solutions[i].seq == outseq)
               ind = i ;
         if (ind < 0)
            break ;
         pending_solutions[ind].report() ;
         pending_solutions.erase(pending_solutions.begin()+ind) ;
         outseq++ ;
      }
   } else {
      if (allsols) {
         int cnt = (sol.length < 1 ? 1 : sol.moves.size() / sol.length) ;
         cout << "Saving " << cnt << " solutions " << sol.seq << endl ;
      } else {
         cout << "Saving solution " << sol.seq << endl ;
      }
      pending_solutions.push_back(sol) ;
   }
}
int parity(const cubepos &cp) {
   int seen = 0 ;
   int rv = 0 ;
   for (int i=0; i<8; i++) {
      if (((seen >> i) & 1) == 0) {
         int cl = 0 ;
         int at = i ;
         while (((seen >> at) & 1) == 0) {
            seen |= 1 << at ;
            cl++ ;
            at = cubepos::corner_perm(cp.c[at]) ;
         }
         if ((cl & 1) == 0)
            rv++ ;
      }
   }
   return rv & 1 ;
}

@ Symmetry code.
This is some stuff for symmetry.

@(nxopt.cpp@>=
int calcsymm(cubepos cp) {
   cubepos v[NMOVES] ;
   cubepos cp2, cp3 ;
   int n = 0 ;
   int r = 0 ;
   int skipped = 0 ;
   for (int f=0; f<FACES; f++) {
      int oldn = n ;
      for (int t=0; t<TWISTS; t++) {
         cp2 = cp ;
         cp2.move(f*TWISTS+t) ;
         cp2.canon_into96(cp3) ;
         int add = 1 ;
         for (int i=0; i<n; i++)
            if (v[i] == cp3) {
               add = 0 ;
               break ;
            }
         if (add)
            v[n++] = cp3 ;
      }
      if (n == oldn) {
         r |= 1<<f ;
         skipped++ ;
      }
   }
   get_global_lock() ;
   cout << "Skipping " << skipped << " faces." << endl ;
   release_global_lock() ;
   return r ;
}

@ More code.

@(nxopt.cpp@>=
struct worker {
   double start ;
   double durr() {
      double now = walltime() ;
      return now - start ;
   }
   void starttimer() {
      start = walltime() ;
   }
   long long probes, evals ;
   int lookup(const cubepos &cp) {
      probes++ ;
      int cori, cperm ;
      getcornercoords(cp, cperm, cori) ;
      firstlev &fl = firstlevs[cori] ;
      int m = fl.v[cperm].m ;
      int ep = getedgecoord(cp, m) ;
      int off = MEMOFFSET*fl.v[cperm].ind+(ep>>4) ;
      int r = (fl.p[off]>>(2*(ep&15))) & 3 ;
      if (r == 0) {
         r = fl.p[off & ~3] & 15 ;
#ifdef TRACE
         cout << " *" << r ;
#endif
         return r ;
      } else {
#ifdef TRACE
         cout << " " << (r + base) ;
#endif
         return r + base ;
      }
   }
   int lookup(const cubepos &cp, int m) {
      probes++ ;
      int cori, cperm ;
      getcornercoords(cp, cperm, cori, m) ;
      firstlev &fl = firstlevs[cori] ;
      int m2 = fl.v[cperm].m ;
      int ep = getedgecoord(cp, cubepos::mm[m][m2]) ;
      int off = MEMOFFSET*fl.v[cperm].ind+(ep>>4) ;
      int r = (fl.p[off]>>(2*(ep&15))) & 3 ;
      if (r == 0) {
         r = fl.p[off & ~3] & 15 ;
#ifdef TRACE
         cout << " *" << r ;
#endif
         return r ;
      } else {
#ifdef TRACE
         cout << " " << (r + base) ;
#endif
         return r + base ;
      }
   }
   int lookup3(const cubepos &cp, int over, int &vals, int skipat, int skipval) {
      int r = 0 ;
      int v1 = (skipat == 0 ? skipval : lookup(cp)) ;
      if (v1 > over)
         return v1 ;
      vals |= v1 ;
      if (v1 == over)
         r |= 256*9 ;
      int v2 = skipval ;
      if (skipat != 1) {
         v2 = lookup(cp, 32) ;
         if (v2 > over)
            return v2 ;
      }
      vals |= v2 << 4 ;
      if (v2 == over)
         r |= 512*9 ;
      if (v1 > v2)
         swap(v1, v2) ;
      int v3 = skipval ;
      if (skipat != 2) {
         v3 = lookup(cp, 16) ;
         if (v3 > over)
            return v3 ;
      }
      vals |= v3 << 8 ;
      if (v3 == over)
         r |= 1024*9 ;
      if (v1 == v2 && v1 == v3 && v1 != 0)
         return r+v1+1 ;
      if (v3 > v2)
         return r+v3 ;
      else
         return r+v2 ;
   }
   int lookup6(const cubepos &cp, int over, int &vals, int skipat, int skipval) {
      evals++ ;
#ifdef TRACE
 cout << " {" << over << "," << skipat << "," << skipval << "}" ;
#endif
      vals = 0 ;
      int rv1 = lookup3(cp, over, vals, skipat, skipval) ;
      int v1 = rv1 & 255 ;
      if (v1 > over)
         return v1 ;
      int vals0 = vals ;
      vals = 0 ;
      cubepos cp2 ;
      cp.invert_into(cp2) ;
      int rv2 = lookup3(cp2, over, vals, skipat-3, skipval) ;
      vals = (vals << 12) | vals0 ;
      int v2 = rv2 & 255 ;
      int r = ((rv1 & -256) << 6) + (rv2 & -256) ;
      if (v1 > v2)
         return r+v1 ;
      else
         return r+v2 ;
   }
   void show(const cubepos &cp, const char *lab, int m) {
      cubepos cp2, cp3, cp4, cp5, cp6 ;
      cp.remap_into(16, cp2) ;
      cp.remap_into(32, cp3) ;
      cp.invert_into(cp4) ;
      cp4.remap_into(16, cp5) ;
      cp4.remap_into(32, cp6) ;
      int tmp ;
      cout << lab << " " << m << " " << lookup(cp) << " " << lookup(cp2) << " "
           << lookup(cp3) << " " << lookup(cp4) << " " << lookup(cp5) << " "
           << lookup(cp6) << " " << lookup6(cp, 100, tmp, -1, 0) << endl ;
   }
   vector<int> front ;
   vector<int> back ;
   cubepos workingcp ;
   solution sol ;
   vector<pair<ll, int> > heads ;
   void topn(const cubepos &cp, ll lfmask, int d) {
#if defined(SLICE) || defined(AXIAL)
#define TOPSIZE 3
#else
#define TOPSIZE 4
#endif
      if (d <= TOPSIZE+base || allsols) {
         recur(cp, lfmask, CANONSEQSTART, (1LL << NMOVES) - 1,
               CANONSEQSTART, d, -1, 0) ;
         return ;
      }
      sort(heads.begin(), heads.end()) ;
      if (heads.size() == 0) {
         for (int m1=0; m1<NMOVES; m1++) {
            if (!((lfmask >> m1) & 1))
               continue ;
            int cs2 = cubepos::next_cs(CANONSEQSTART, m1) ;
            ll lfmask2 = cubepos::cs_mask(cs2) ;
            for (int m2=0; m2<NMOVES; m2++) {
               if (!((lfmask2 >> m2) & 1))
                  continue ;
               int cs3 = cubepos::next_cs(CANONSEQSTART, m2) ;
               ll lfmask3 = cubepos::cs_mask(cs3) ;
               for (int m3=0; m3<NMOVES; m3++) {
                  if (!((lfmask3 >> m3) & 1))
                     continue ;
#if !defined(SLICE) && !defined(AXIAL)
                  int cs4 = cubepos::next_cs(CANONSEQSTART, m3) ;
                  ll lfmask4 = cubepos::cs_mask(cs4) ;
                  for (int m4=0; m4<NMOVES; m4++) {
                     if (!((lfmask4 >> m4) & 1))
                        continue ;
                     heads.push_back(make_pair(0LL,
                            m1 + NMOVES * (m2 + NMOVES * (m3 + NMOVES * m4)))) ;
#else
                     heads.push_back(make_pair(0LL,
                            m1 + NMOVES * (m2 + NMOVES * m3))) ;
#endif
#if !defined(SLICE) && !defined(AXIAL)
                  }
#endif
               }
            }
         }
      }
      for (int ohi=heads.size()-1; ohi>=0; ohi--) {
         front.clear() ;
         back.clear() ;
         ll oevals = evals ;
         int ms = heads[ohi].second ;
         cubepos cp2 = cp ;
         int cs = CANONSEQSTART ;
         int mv = ms % NMOVES ;
#if !defined(SLICE) && !defined(AXIAL)
         front.push_back(mv) ;
         cp2.movepc(mv) ;
         cs = cubepos::next_cs(cs, mv) ;
         ms /= NMOVES ;
         mv = ms % NMOVES ;
#endif
         front.push_back(mv) ;
         cp2.movepc(mv) ;
         cs = cubepos::next_cs(cs, mv) ;
         ms /= NMOVES ;
         mv = ms % NMOVES ;
         front.push_back(mv) ;
         cp2.movepc(mv) ;
         cs = cubepos::next_cs(cs, mv) ;
         ms /= NMOVES ;
         mv = ms % NMOVES ;
         front.push_back(mv) ;
         cp2.movepc(mv) ;
         cs = cubepos::next_cs(cs, mv) ;
         int tt = recur(cp2, cubepos::cs_mask(cs), cs, (1LL << NMOVES) - 1,
                        CANONSEQSTART, d-TOPSIZE, -1, 0) ;
         if (tt == 0 && !allsols)
            return ;
         oevals = evals - oevals ;
         heads[ohi].first = oevals / order_mult[cs] ;
      }
   }
   int recur(const cubepos &cp, long long lfmask, int fprev, long long rlfmask,
             int rprev, int togo, int skipat, int skipval) {
      if (togo == 0) {
#ifdef TRACE
         cout << endl << "Possible solution" ;
 for (int i=0; i<8; i++) cout << " " << (int)(cp.c[i]) ;
 for (int i=0; i<12; i++) cout << " " << (int)(cp.e[i]) ;
 cout << endl ;
#endif
         if (cp == identity_cube) {
            testit(workingcp, sol) ;
            return 0 ;
         } else {
            return 1 ;
         }
      }
//    show(cp, "togo", togo) ;
      int vals = 0 ;
      int lr = lookup6(cp, togo, vals, skipat, skipval) ;
      int v = lr & 255 ;
      if (v > togo)
         return v ;
      long long rfmask = expandm[(lr >> 14) & 07] & rlfmask ;
      long long fmask = 0 ;
      cubepos cp2 ;
      togo-- ;
      int dir = 0 ;
      if (bidir) {
         fmask = expandm[(lr >> 8) & 07] & lfmask ;
         dir = popcount64(rfmask) - popcount64(fmask) ;
         if (dir == 0) {
            dir = (vals & 15) + ((vals >> 4) & 15) + ((vals >> 8) & 15)
             - ((vals >> 12) & 15) - ((vals >> 16) & 15) - ((vals >> 20) & 15) ;
         }
      }
      if (dir > 0) {
         for (int mv=0; mv<NMOVES; mv++) {
            if ((fmask >> mv) == 0)
               break ;
            if (0 == ((fmask >> mv) & 1))
               continue ;
            cp2 = cp ;
            cp2.movepc(mv) ;
            skipat = 3 + skipata[mv] ;
            skipval = (vals >> (4 * skipat)) & 15 ;
            int ncs = cubepos::next_cs(fprev, mv) ;
            front.push_back(mv) ;
            int tt = recur(cp2, cubepos::cs_mask(ncs),
                           ncs, rlfmask, rprev, togo, skipat, skipval) ;
            if (tt == 0 && !allsols)
               return 0 ;
            front.pop_back() ;
#if defined(HALF) || defined(SLICE)
            if (tt > togo+1)
               fmask &= ~(7LL << (mv / TWISTS * TWISTS)) ;
#endif
#ifdef AXIAL
            if (tt > togo+1)
               fmask &= ~(07007007007007LL << (3*skipata[mv])) ;
#endif
         }
      } else {
         for (int mv=0; mv<NMOVES; mv++) {
            if ((rfmask >> mv) == 0)
               break ;
            if (0 == ((rfmask >> mv) & 1))
               continue ;
            cp2 = cp ;
            cp2.move(mv) ;
            skipat = skipata[mv] ;
            skipval = (vals >> (4 * skipat)) & 15 ;
            int ncs = cubepos::next_cs(rprev, mv) ;
            back.push_back(mv) ;
            int tt = recur(cp2, lfmask, fprev,
                           cubepos::cs_mask(ncs), ncs, togo, skipat, skipval) ;
            if (tt == 0 && !allsols)
               return 0 ;
            back.pop_back() ;
#if defined(HALF) || defined(SLICE)
            if (tt > togo+1)
               rfmask &= ~(7LL << (mv / TWISTS * TWISTS)) ;
#endif
#ifdef AXIAL
            if (tt > togo+1)
               rfmask &= ~(07007007007007LL << (3*skipata[mv])) ;
#endif
         }
      }
      if (v == 0)
         return 1 ;
      else
         return v ;
   }
   int solve(const cubepos &cp, int seq) {
      probes = 0 ;
      evals = 0 ;
      front.clear() ;
      back.clear() ;
      sol.length = -1 ;
      int tmp = 0 ;
#ifdef QUARTER
      int cparity = parity(cp) ;
#endif
      int d = lookup6(cp, 100, tmp, -1, 0) ;
      if (d < mindepth)
         d = mindepth ;
      long long lfmask = (1LL << NMOVES) - 1 ;
      if (symmetry)
         lfmask = calcsymm(cp) ;
      heads.clear() ;
      for (int d=lookup6(cp, 100, tmp, -1, 0); d <= maxdepth; d++) {
#ifdef QUARTER
         if ((d ^ cparity) & 1)
            continue ;
#endif
#ifndef QUARTER
         if (leftover)
            lfmask &= expandm[1] ;
#endif
         topn(cp, lfmask, d) ;
         get_global_lock() ;
         cout << "Finished searching depth " << d << " pos " << seq
             << " in " << probes << " probes " << evals << " evals." << endl ;
         release_global_lock() ;
         if (sol.length >= 0)
            return d ;
      }
      heads.clear() ;
      return -1 ;
   }

@ More code.

@(nxopt.cpp@>=
   void testit(cubepos cp, solution &sol) {
      int len = 0 ;
      for (int i=0; i<(int)front.size(); i++) {
        int mv = front[i] ;
        cp.movepc(mv) ;
        sol.moves.push_back(mv) ;
        len++ ;
      }
      for (int i=(int)back.size()-1; i>=0; i--) {
        int mv = cubepos::inv_move[back[i]] ;
        cp.movepc(mv) ;
        sol.moves.push_back(mv) ;
        len++ ;
      }
      if (!(cp == identity_cube))
         error("! solver failed") ;
      sol.length = len ;
      savestats(sol) ;
   }
   void savestats(solution &sol) {
      sol.evals = evals ;
      sol.probes = probes ;
      sol.duration = durr() ;
   }
   void dowork() {
#ifdef TESTTABLE
      {
         ull cnt[30], cnt6[30] ;
         for (int i=0; i<30; i++)
            cnt[i] = cnt6[i] = 0 ;
         cubepos cp ;
         for (int i=0; i<1000; i++)
            cp.move(random_move()) ;
         const int TESTS = 100000000 ;
         duration() ;
         int tmp ;
         for (int i=0; i<TESTS; i++) {
            cnt[lookup(cp)]++ ;
            probes-- ;
            cnt6[lookup6(cp, base+3, tmp, -1, 0)&255]++ ;
            cp.move(random_move()) ;
         }
         cout << "Did " << TESTS << " probes " << probes << " in " << duration() << endl ;
         double sum = 0, sum6 = 0 ;
         for (int i=0; i<30; i++)
            if (cnt[i] != 0 || cnt6[i] != 0) {
               cout << i << " " << cnt[i] << " " << cnt6[i] << endl ;
               sum += i * cnt[i] ;
               sum6 += i * cnt6[i] ;
            }
         cout << "Average " << (sum/TESTS) << " " << (sum6/TESTS) << endl ;
         exit(0) ;
      }
#endif
      while (1) {
         solution soli ;
         cubepos cp ;
         int gotwork = 0 ;
         get_global_lock() ;
         gotwork = getwork(cp) ;
         release_global_lock() ;
         if (!gotwork)
            return ;
         sol = soli ;
         sol.seq = gotwork ;
         starttimer() ;
         workingcp = cp ;
         sol.length = solve(cp, gotwork) ;
         savestats(sol) ;
         get_global_lock() ;
         report(sol) ;
         release_global_lock() ;
      }
   }
   int buf[64] ; // padding so workers don't share cache lines
} workers[MAX_THREADS] ;
void parseposition(cubepos &cp, char *s) {
  cp = identity_cube ;
  if (cp.parse_Singmaster(s) == 0)
    return ;
  cp = identity_cube ;
  const char *faces = "UuFfRrDdBbLl" ;
  while (*s) {
    while (*s && *s <= ' ')
       s++ ;
    if (*s == 0)
       break ;
    const char *fp = index(faces, *s) ;
    if (fp == 0) {
       cerr << "Illegal face in " << *s << endl ;
       exit(10) ;
    }
    int f = (fp - faces) / 2 ;
    s++ ;
    switch (*s) {
case '-':
case '3':
       cp.movepc(f*TWISTS) ;
case '2':
       cp.movepc(f*TWISTS) ;
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
char myline[1001] ;
int fromstdin = 0 ;
cubepos seqcp ;
int maxrandom = 2000000000 ;
double maxsecs = 1e20 ;
double starttime ;
int inputseq = 1 ;
int getwork(cubepos &cp) {
   if (walltime() > starttime + maxsecs)
      exit(0) ;
   if (fromstdin) {
      if (fgets(myline, 1000, stdin) == 0)
         return 0 ;
      parseposition(cp, myline) ;
   } else {
      if (maxrandom-- <= 0)
         return 0 ;
      cp = seqcp ;
      int mv = (int)(NMOVES*erand48(entropy)) ;
      seqcp.movepc(mv) ;
   }
   return inputseq++ ;
}
void *pworker(void *s) {
   seed48(entropy) ; // needed only for cygwin
   worker *w = (worker*)s ;
   w->dowork() ;
   return 0 ;
}
int main(int argc, char *argv[]) {
   cout << "This is " << argv[0] << " metric " << metric << endl ;
   for (int i=0; i<argc; i++)
      cout << " " << argv[i] ;
   cout << endl << flush ;
   base = BASE ;
   srand48(0) ; // needed only for cygwin
   seqcp = identity_cube ;
   init() ;
   while (argc > 1 && argv[1][0] == '-' && argv[1][1]) {
      argc-- ;
      argv++ ;
      switch (argv[0][1]) {
case 'b':
         sscanf(argv[1], "%d", &base) ;
         argc-- ;
         argv++ ;
         break ;
case 'c':
         sscanf(argv[1], "%d", &maxrandom) ;
         argc-- ;
         argv++ ;
         break ;
case 'T':
         maxsecs = atof(argv[1]) ;
         argc-- ;
         argv++ ;
         break ;
#ifdef THREADS
case 't':
         sscanf(argv[1], "%d", &numthreads) ;
         argc-- ;
         argv++ ;
         if (numthreads < 0)
            numthreads = 1 ;
         if (numthreads > MAX_THREADS)
            numthreads = MAX_THREADS ;
         break ;
#endif
case 'n':
         mindepth = atol(argv[1]) ;
         argc-- ;
         argv++ ;
         break ;
case 'x':
         maxdepth = atol(argv[1]) ;
         argc-- ;
         argv++ ;
         break ;
case 's':
         symmetry++ ;
         break ;
#ifndef QUARTER
case 'l':
         leftover++ ;
         break ;
#endif
case 'U':
         bidir = 0 ;
         break ;
case 'A':
         allsols++ ;
         break ;
default:
         error("! bad arg") ;
      }
   }
   if (allsols && bidir)
      cout << "Warning: using bidirectional and all simultaneously may give "
           << "unexpected results." << endl ;
   // calculate order multipliers for canonical sequences.  Proportional
   // to the number of sequences starting in a given state.
   for (int i=0; i<CANONSEQSTATES; i++)
      order_mult[i] = 1 ;
   double norder_mult[CANONSEQSTATES] ;
   for (int i=0; i<20; i++) {
      for (int cs=0; cs<CANONSEQSTATES; cs++)
         norder_mult[cs] = 0 ;
      for (int cs=0; cs<CANONSEQSTATES; cs++) {
         ll mvmask = cubepos::cs_mask(cs) ;
         for (int mv=0; mv<NMOVES; mv++) {
            if (!((mvmask >> mv) & 1))
               continue ;
            norder_mult[cs] += order_mult[cubepos::next_cs(cs, mv)] ;
         }
      }
      double sum = 0 ;
      for (int cs=0; cs<CANONSEQSTATES; cs++)
         sum += norder_mult[cs] ;
      for (int cs=0; cs<CANONSEQSTATES; cs++)
         order_mult[cs] = norder_mult[cs] / sum ;
   }
   starttime = walltime() ;
   if (!readtab())
      generatetab() ;
   seed48(entropy) ;
#ifdef TESTTABLE
   workers[0].dowork() ;
#endif
   fromstdin = (argc > 1 && strcmp(argv[1], "-") == 0) ;
#ifdef THREADS
   pthread_t p_thread[MAX_THREADS] ;
   for (int ti=0; ti<numthreads; ti++) {
      pthread_create(&(p_thread[ti]), NULL, pworker, workers+ti) ;
   }
   for (int ti=0; ti<numthreads; ti++) {
      pthread_join(p_thread[ti], 0) ;
   }
#else
   workers[0].dowork() ;
#endif
}
