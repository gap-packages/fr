#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <mt.h>

#define LEVEL 6
#define FIRST_GRIGORCHUK_GROUP
#define ERSCHLER_GROUP_
#define GRIGORCHUK_OVERGROUP_

const unsigned long left[] =
  { 0xffff0000, 0xff00ff00, 0xf0f0f0f0, 0xcccccccc, 0xaaaaaaaa, 0x00000000 };
const unsigned long right[] =
  { 0x0000ffff, 0x00ff00ff, 0x0f0f0f0f, 0x33333333, 0x55555555, 0xffffffff };
const int shift[] =
  { 16, 8, 4, 2, 1, 0 };
const int xshift[] =
  { 16, 8, 4, 2, 1, 1 };

unsigned long swap (unsigned long w, unsigned long mask, int level) {
  return (w & ~mask) | ((w & mask & left[level]) >> shift[level]) | ((w & mask & right[level]) << shift[level]) ;
}

#if LEVEL == 6
typedef unsigned long long htype;
#else
typedef unsigned long htype;
#endif

struct treeaut {
  unsigned long portrait[LEVEL];

  treeaut () {}
  treeaut (char t) {
    switch (t) {
    case '1':
      for (int i = 0; i < LEVEL; i++) portrait[i] = 0;
      break;
    case 'a':
      portrait[0] = ~0;
      for (int i = 1; i < LEVEL; i++) portrait[i] = 0;
      break;
#ifdef FIRST_GRIGORCHUK_GROUP
    case 'b': case 'c': case 'd':
      {
	portrait[0] = 0;
	int mask = ~0;
	for (int i = 1; i < LEVEL; i++) {
	  if ((i-1 + t -'d') % 3)
	    portrait[i] = mask & left[i-1];
	  else
	    portrait[i] = 0;
	  mask &= right[i-1];
	}
	return;
      }      
#elif defined(ERSCHLER_GROUP)
    case 'b': case 'c':
      {
	portrait[0] = 0;
	int mask = ~0;
	for (int i = 1; i < LEVEL; i++) {
	  if ((i-1 + t-'c') % 2)
	    portrait[i] = mask & left[i-1];
	  else
	    portrait[i] = 0;
	  mask &= right[i-1];
	}
	return;
      }
    case 'd':
      {
	portrait[0] = 0;
	int mask = ~0;
	for (int i = 1; i < LEVEL; i++) {
	  portrait[i] = mask & left[i-1];
	  mask &= right[i-1];
	}
	return;
      }
#elif defined(GRIGORCHUK_OVERGROUP)
    case 'b': case 'c': case 'd':
    case 'e': case 'f': case 'g':
    case 'h':
      {
	portrait[0] = 0;
	int mask = ~0;
	for (int i = 1; i < LEVEL; i++) {
	  if ((1 << ((i-1)%3)) & (t-'b'+1))
	    portrait[i] = mask & left[i-1];
	  else
	    portrait[i] = 0;
	  mask &= right[i-1];
	}
	return;
      }      
#else
#error No group defined
#endif
    }
  }
  int activity() { return portrait[0] == 0; }
  bool trivial() { for (int i = 0; i < LEVEL; i++) if (portrait[i]) return false; return true; }
  htype packed (void) {
    unsigned long w = 0;
    for (int i = 0; i < LEVEL && i <= 4; i++)
#ifdef GRIGORCHUK_OVERGROUP
      {
	for (unsigned long j = 0x80000000; j; j >>= shift[i], j >>= xshift[i]) {
	  if (i == 4 && (j == 0x80000000 || j == 0x00800000 || j == 0x00008000))
	    continue;
	  w <<= 1;
	  if (j & portrait[i])
	    w |= 1;
	}
      }
#else
      w ^= portrait[i] & right[i];
#endif
    if (LEVEL == 6)
      return ((htype) w) << 32 | portrait[5];
    else
      return w;
  }
  treeaut &operator *= (treeaut x) {
    for (int i = 0; i < LEVEL; i++) {
      for (int j = i; j < LEVEL; j++)
	x.portrait[j] = swap(x.portrait[j], portrait[i], i);
      portrait[i] ^= x.portrait[i];
    }
    return *this;
  }
  treeaut operator * (treeaut &x) {
    treeaut n;

    n = *this;
    n *= x;
    return n;
  }
  treeaut operator -(void) {
    treeaut n = *this;
    for (int i = 0; i < LEVEL; i++) {
      for (int j = 0; j < LEVEL; j++)
	n.portrait[j] = swap(n.portrait[j], n.portrait[i], i);
    }
    return n;
  }
};

std::ostream &operator <<(std::ostream &s, treeaut t)
{
  s << "[";
  for (int i = 0; i < LEVEL; i++) {
    for (unsigned long j = 0x80000000; j ; j >>= shift[i], j >>= xshift[i]) {
      if (j != 0x80000000 && i <= 4)
	s << ' ';
      if (j & t.portrait[i]) s << '@'; else s << '.';
    }
    if (i == LEVEL-1) s << ']'; else s << '|';
  }
  return s;
}

#define HSIZE 0x02000000
#if LEVEL == 6
#define HFREE 0xffffffffffffffffLL
#else
#define HFREE 0xffffffff
#endif

#ifndef PERFECT_HASHING
htype hkey[HSIZE];
#endif
char hdata[HSIZE];

htype hashkey (treeaut &t)
{
#ifdef PERFECT_HASHING
  return t.packed();
#else
  htype p = t.packed();
  unsigned long x = p & (HSIZE-1);

  while (hkey[x] != p) {
    if (hkey[x] == HFREE) { hkey[x] = p; return x; }
    x = (x + 1) & (HSIZE-1);
  }
  return x;
#endif
}

void print_key (unsigned long k) {
  std::cout << "[";
  for (int j = 0; j < 28; j++, k <<= 1) {
    if (k & 0x08000000) std::cout << "@"; else std::cout << ".";
    switch (j) {
    case 0: case 2: case 6:
      std::cout << "|"; break;
    case 27: break;
    case 14:
      std::cout << "|? "; break;
    case 17: case 20:
      std::cout << " ?";
    default: std::cout << " ";
    }
  }
  std::cout << "]";
}

treeaut trivial('1'), A('a'), B('b'), C('c'), D('d');
#if defined(GRIGORCHUK_OVERGROUP)
treeaut E('e'), F('f'), G('g'), H('h');
#endif

int level, maxdepth, diameter;

void fillin (treeaut t, int depth)
{
  //  for (int i = 0; i < depth; i++) std::cout << "  "; cout << depth << t << "\n";
  unsigned long k = hashkey(t);
  if (hdata[k] == -1 || hdata[k] > depth) {
    hdata[k] = depth;
    // std::cout << t << " depth: " << depth << "\n";
    if (depth == maxdepth) return;
    fillin (t * A, depth+1);
    fillin (t * B, depth+1);
    fillin (t * C, depth+1);
    fillin (t * D, depth+1);
#if defined(GRIGORCHUK_OVERGROUP)
    fillin (t * E, depth+1);
    fillin (t * F, depth+1);
    fillin (t * G, depth+1);
    fillin (t * H, depth+1);
#endif
  }
}

void printgeod (treeaut t, int depth)
{
  static char mem[256];

  unsigned long k = hashkey(t);
  if (hdata[k] != depth) return;

  if (depth == maxdepth) {
    mem[depth] = 0;
    std::cout << mem+level << " " << t << "\n";
    // std::cout << "++++++++ "; print_key (k); std::cout << "\n";
    return;
  }
  mem[depth] = 'a'; printgeod (t * A, depth+1);
  mem[depth] = 'b'; printgeod (t * B, depth+1);
  mem[depth] = 'c'; printgeod (t * C, depth+1);
  mem[depth] = 'd'; printgeod (t * D, depth+1);
#if defined(GRIGORCHUK_OVERGROUP)
  mem[depth] = 'e'; printgeod (t * E, depth+1);
  mem[depth] = 'f'; printgeod (t * F, depth+1);
  mem[depth] = 'g'; printgeod (t * G, depth+1);
  mem[depth] = 'h'; printgeod (t * H, depth+1);
#endif
}

int main (int argc, char *argv[]) {
  memset (hdata, -1, sizeof hdata);
#ifndef PERFECT_HASHING
  memset (hkey, -1, sizeof hkey);
#endif

  if (argc != 2) {
    std::cerr << "Use: " << argv[0] << " <maxdepth>\n";
    return -1;
  }
  maxdepth = atoi(argv[1]);

#if FAST_SCAN_BY_5
  for (int i = maxdepth - (maxdepth % 5); i >= 0; i -= 5)
    fillin (trivial, i);
#else
  for (level = maxdepth; level >= 0; level --) {
    fillin (trivial, level);
    
    int hist[256], missed = 0x400000;
    memset (hist, 0, sizeof hist);

    for (int i = 0; i < HSIZE; i++) {
      if (hdata[i] >= 0)
	hist[hdata[i]]++, missed--;
    }
    
    for (int i = level; i < 256 && hist[i]; i++)
      std::cout << "Len " << i-level << ": " << hist[i] << "\n";
    std::cout << "Missed: " << missed << "\n";
    std::cout.flush();
    if (missed == 0)
      break;
  }
#endif
  int hist[256];
  memset (hist, 0, sizeof hist);

  for (int i = 0; i < HSIZE; i++) {
    if (hdata[i] >= 0)
      hist[hdata[i]]++;
#ifdef PRINT_KEYS
    if (hdata[i] == 8) {
      std::cout << "???????? "; print_key(i); std::cout << "\n";
    }
#endif
#ifdef PRINT_KEYS
    if (hdata[i] == maxdepth) {
      if (hkey[i] & 0xaaaaaaaa)
	std::cout << "KEY " << std::hex << hkey[i] << std::dec << "\n";
      else {
	std::cout << "KEY " << std::hex << hkey[i] << std::dec << "CODE ";
	for (int j = 0x40000000; j; j >>= 2)
	  if (hkey[i] & j) std::cout << "*"; else std::cout << " ";
	std::cout << "\n";
      }
    }
#endif
  }

  for (int i = level; i < 256 && hist[i]; i++)
    std::cout << "Len " << i-level << ": " << hist[i] << "\n";
  std::cout.flush();
  
  printgeod (trivial, level);

  return 0;
}
/*
 *LEVEL 1: [1,1]
 *LEVEL 2: [1,2,2,2,1]
 *LEVEL 3: [1,4,6,12,16,20,24,28,17]
 *LEVEL 4: [1,4,6,12,17,28,40,68,92,132,
 *          168,236,291,356,436,556,493,308,282,268,
 *	    182,60,35,20,5]
 *LEVEL 5: [1,4,6,12,17,28,40,68,95,156,
 *          216,356,488,772,1054,1660,2218,3332,4450,6700,
 *	    8773,12716,16538,23924,30161,40820,51444,69836,85796,112372,
 *	    138708,182684,202815,216828,240300,267940,278748,282052,290753,293308,
 *	    271007,228380,210627,186868,148625,95804,73495,53044,31765,13004,
 *	    7305,4020,1495,368,193,96,19]
 *LEVEL 6: [1,4,6,12,17,28,40,68,95,156,
 *          216,356,488,772,1054,1660,2240,3448,4642,7128,
 *          9518,14392,19186,28984,38237,56452,74422,109676,143511,209636,
 *          274674,400028,...518950...]
 */
