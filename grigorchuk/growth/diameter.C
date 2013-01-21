#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include "lz4-read-only/lz4.h"

#ifdef TEST
#define BUFSIZE 0x4000L
#define LEVEL 5
#define GROUPSIZE (1L << 22)
#define MAXTREE (1L << 18)
const char prefix[] = "/tmp/datatest_%03d";
#define FULLGROUP
#else
#define BUFSIZE 0x40000000L
#define LEVEL 6
#define GROUPSIZE (1L << 42)
#define MAXTREE (1L << 30)
const char prefix[] = "/bigwork/ncablaur/data6_%03d";
#define FULLGROUP
#endif

#define FIRST_GRIGORCHUK_GROUP
//#define ERSCHLER_GROUP
//#define GRIGORCHUK_OVERGROUP

#if LEVEL <= 6
typedef unsigned long activitymask; // in portrait
typedef unsigned long htype; // (perfect) hash key, needs double precision
const activitymask right[] =
  { 0xffff0000, 0xff00ff00, 0xf0f0f0f0, 0xcccccccc, 0xaaaaaaaa, 0x00000000 };
const activitymask left[] =
  { 0x0000ffff, 0x00ff00ff, 0x0f0f0f0f, 0x33333333, 0x55555555, 0xffffffff };
const activitymask copy[] =
  { 0xffffffff, 0x0000ffff, 0x000000ff, 0x0000000f, 0x00000003, 0x00000001 };
const activitymask topmask =
  0x80000000;
const int halfspread[] =
  { 16, 8, 4, 2, 1, 0 };
const int spread[] =
  { 32, 16, 8, 4, 2, 1 };
#else
typedef unsigned long activitymask;
typedef __int128_t htype;
const activitymask right[] =
  { 0xffffffff00000000L, 0xffff0000ffff0000L, 0xff00ff00ff00ff00L, 0xf0f0f0f0f0f0f0f0L, 0xccccccccccccccccL, 0xaaaaaaaaaaaaaaaaL, 0x0000000000000000L };
const activitymask left[] =
  { 0x00000000ffffffffL, 0x0000ffff0000ffffL, 0x00ff00ff00ff00ffL, 0x0f0f0f0f0f0f0f0fL, 0x3333333333333333L, 0x5555555555555555L, 0xffffffffffffffffL };
const activitymask copy[] =
  { 0xffffffffffffffffL, 0x00000000ffffffffL, 0x000000000000ffffL, 0x00000000000000ffL, 0x000000000000000fL, 0x0000000000000003L, 0x0000000000000001L };
const activitymask topmask =
  0x8000000000000000L;
const int halfspread[] =
  { 32, 16, 8, 4, 2, 1, 0 };
const int spread[] =
  { 64, 32, 16, 8, 4, 2, 1 };
#endif

const int pow2[] =
  { 1, 2, 4, 8, 16, 32, 64, 128 };

#pragma omp threadprivate(right,left,copy,topmask,halfspread,spread,pow2)

// in w, swap bits selected by mask between positions left[level] and right[level] 
activitymask swap (activitymask w, activitymask mask, int level) {
  return (w & ~mask) | ((w & mask & left[level]) << halfspread[level]) | ((w & mask & right[level]) >> halfspread[level]) ;
}

struct treeaut {
  activitymask portrait[LEVEL];
  // we store an automorphism as a sequence of bitsets, one for each level.
  // level 0 corresponds to the top of the tree.
  // in portrait[l], all bits == i mod 2^l are equal, and give the activity
  // on vertex i at level l.
  int activity() { return portrait[0] == 0; }
  bool trivial() { for (int i = 0; i < LEVEL; i++) if (portrait[i]) return false; return true; }
  htype pack (void) const {
    htype key = 0;
    for (int i = LEVEL-1; i >= 0; i--) {
      activitymask s = topmask;
      for (int j = pow2[i]-1; j >= 0; j--, s >>= spread[i]) {
#ifdef FIRST_GRIGORCHUK_GROUP
	if (i >= 3 && (j%8 == 3 || j%8 == 5 || j%8 == 7))
	  continue;
	key <<= 1;
	if (portrait[i] & s)
	  key |= 1;
#else
	std::cerr << "No packer!\n";
#endif
      }
    }
    return key;
  }
  void unpack (htype key) {
    for (int i = 0; i < LEVEL; i++) {
      portrait[i] = 0;
      activitymask s = copy[i];
      for (int j = 0; j < pow2[i]; j++, s <<= spread[i]) {
#ifdef FIRST_GRIGORCHUK_GROUP
	if (i >= 3 && (j%8 == 3 || j%8 == 5 || j%8 == 7)) {
	  if (j%8 == 3)
	    if (((portrait[i] << spread[i]) ^ (portrait[i] << (2*spread[i])) ^ (portrait[i] << (3*spread[i])) ^ (portrait[i-1]) ^ (portrait[i-1] << spread[i-1]) ^ (portrait[i-2] >> spread[i-2])) & s)
	      portrait[i] |= s;
	  if (j%8 == 5)
	    if (((portrait[i] << spread[i]) ^ (portrait[i] << (4*spread[i])) ^ (portrait[i] << (5*spread[i])) ^ (portrait[i-1] >> spread[i-1]) ^ (portrait[i-1] << spread[i-1]) ^ (portrait[i-2] & (portrait[i-2] << spread[i-2]))) & s)
	      portrait[i] |= s;
	  if (j%8 == 7)
	    if (((portrait[i] << spread[i]) ^ (portrait[i] << (2*spread[i])) ^ (portrait[i] << (3*spread[i])) ^ (portrait[i-1]) ^ (portrait[i-1] << spread[i-1]) ^ (portrait[i-2] << spread[i-2])) & s)
	      portrait[i] |= s;
	} else {
	  if (key & 1)
	    portrait[i] |= s;
	  key >>= 1;
	}
#else
	std::cerr << "No unpacker!\n";
#endif
      }
    }
  }
  void generator (const char type) {
    switch (type) {
    case '1':
      for (int i = 0; i < LEVEL; i++) portrait[i] = 0;
      break;
    case 'a':
      portrait[0] = copy[0];
      for (int i = 1; i < LEVEL; i++) portrait[i] = 0;
      break;
#ifdef FIRST_GRIGORCHUK_GROUP
    case 'b': case 'c': case 'd':
      {
	portrait[0] = 0;
	activitymask mask = copy[0];
	for (int i = 1; i < LEVEL; i++) {
	  if ((i-1 + type -'d') % 3)
	    portrait[i] = mask & left[i-1];
	  else
	    portrait[i] = 0;
	  mask &= right[i-1];
	}
	break;
      }      
#elif defined(ERSCHLER_GROUP)
    case 'b': case 'c': case 'd':
      {
	portrait[0] = 0;
	activitymask mask = copy[0];
	for (int i = 1; i < LEVEL; i++) {
	  if (((i-1 + type-'c') % 2) || (type == 'd'))
	    portrait[i] = mask & left[i-1];
	  else
	    portrait[i] = 0;
	  mask &= right[i-1];
	}
	break;
      }
#elif defined(GRIGORCHUK_OVERGROUP)
    case 'b': case 'c': case 'd':
    case 'e': case 'f': case 'g':
    case 'h':
      {
	portrait[0] = 0;
	activitymask mask = copy[0];
	for (int i = 1; i < LEVEL; i++) {
	  if ((1 << ((i-1)%3)) & (type-'b'+1))
	    portrait[i] = mask & left[i-1];
	  else
	    portrait[i] = 0;
	  mask &= right[i-1];
	}
	break;
      }      
#else
#error No group defined
#endif
    }
  }
};

treeaut operator *(treeaut a, treeaut b) {
  for (int i = 0; i < LEVEL; i++) {
    for (int j = i; j < LEVEL; j++)
      b.portrait[j] = swap(b.portrait[j], a.portrait[i], i);
    a.portrait[i] ^= b.portrait[i];
  }
  return a;
}

treeaut operator -(treeaut a) { // really, inverse
  for (int i = 0; i < LEVEL; i++) {
    for (int j = 0; j < LEVEL; j++)
      a.portrait[j] = swap(a.portrait[j], a.portrait[i], i);
  }
  return a;
}

std::ostream &operator <<(std::ostream &stream, treeaut t)
{
  stream << "[";
  for (int i = 0; i < LEVEL; i++) {
    activitymask s = 1;
    for (int j = 0; j < pow2[i]; j++, s <<= spread[i]) {
      if (j != 0 && i <= 4)
	stream << ' ';
      if (s & t.portrait[i]) stream << '@'; else stream << '.';
    }
    if (i == LEVEL-1) stream << ']'; else stream << '|';
  }
  return stream;
}

#ifdef FULLGROUP

struct booltable {
  unsigned char *data;
  long len;
  unsigned long bitcount[256];

  booltable(long setlen) {
    len = setlen;
    data = new unsigned char[setlen/8];
#pragma omp parallel for
    for (long i = 0; i < setlen/8; i++)
      data[i] = 0;
    for (int i = 0; i < 256; i++) {
      bitcount[i] = 0;
      for (unsigned char c = 0x80; c; c >>= 1)
	if (i & c)
	  bitcount[i]++;
    }
  }

  bool get (unsigned long pos) { return (data[pos/8] & pow2[pos%8]) != 0; }

  void set0 (unsigned long pos) {
#pragma omp atomic
    data[pos/8] &= ~pow2[pos%8];
  }

  void set1 (unsigned long pos) {
#pragma omp atomic
    data[pos/8] |= pow2[pos%8];
  }

  void set (unsigned long pos, bool b) {
    if (b) set1(pos); else set0(pos);
  }

  bool operator [](unsigned long pos) { return get(pos); }

  long size (void) {
    long c = 0;
#pragma omp parallel for reduction(+:c)
    for (long i = 0; i < len/8; i++)
      c += bitcount[data[i]];
    return c;
  }

  void save(char *s) {
    int f = open(s, O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (f < 0) {
      std::cerr << "Error opening file " << s << " for writing\n";
      exit(1);
    }
    char *buf = new char[LZ4_compressBound(BUFSIZE)];
    for (long g = 0; g < len/8; g += BUFSIZE) {
      int len = LZ4_compress((char *) data+g, buf, BUFSIZE);
      if (write (f, &len, 4) != 4) {
	std::cerr << "Error writing data to file " << s << "\n";
	exit(1);
      }
      if (write (f, buf, len) != len) {
	std::cerr << "Error writing data to file " << s << "\n";
	exit(1);
      }
    }
    delete buf;
    close(f);
  }

  void load(char *s) {
    int f = open(s, O_RDONLY);
    if (f < 0) {
      std::cerr << "Error opening file " << s << " for reading\n";
      exit(1);
    }
    char *buf = new char[LZ4_compressBound(BUFSIZE)];
    for (long g = 0; g < len/8; g += BUFSIZE) {
      int len;      
      if (read (f, &len, 4) != 4) {
	std::cerr << "Error reading data from file " << s << "\n";
	exit(1);
      }
      if (read (f, buf, len) != len) {
	std::cerr << "Error reading data from file " << s << "\n";
	exit(1);
      }
      LZ4_uncompress(buf, (char *) data+g, BUFSIZE);
    }
    close(f);
  }
};

typedef std::set<htype> hashset;

treeaut trivial, A, B, C, D;
#pragma omp threadprivate(trivial,A,B,C,D)

#if defined(GRIGORCHUK_OVERGROUP)
treeaut E, F, G, H;
#pragma omp threadprivate(E,F,G,H)
#endif

int main(int argc, char *argv[]) {
  char logname[100];
  sprintf(logname, "log-%d", getpid());
  std::ofstream output(logname, std::ios::out);

#pragma omp parallel
  {
    A.generator('a'), B.generator('b'), C.generator('c'), D.generator('d');
#if defined(GRIGORCHUK_OVERGROUP)
    E.generator('e'), F.generator('f'), G.generator('g'), H.generator('h');
#endif
  }

  hashset table[2];

  long count = 0, oldcount, diam = 0;
  bool resume = (argc > 1);

  if (resume)
    diam = atoi(argv[1]);
    
  if (!resume) {
    table[0].insert(trivial.pack());
    table[1].insert(trivial.pack());

    for (; ; diam++) {
      oldcount = count;
      count = table[diam%2].size();

      std::cout << diam << "\t" << count-oldcount << std::endl;
      output << diam << "\t" << count-oldcount << std::endl;
      output.flush();
    
      if (count > MAXTREE)
	break;

      for (hashset::iterator g = table[diam%2].begin(); g != table[diam%2].end(); g++) {
	treeaut t;
	t.unpack(*g);
	table[(diam+1)%2].insert((t*A).pack());
	table[(diam+1)%2].insert((t*B).pack());
	table[(diam+1)%2].insert((t*C).pack());
	table[(diam+1)%2].insert((t*D).pack());
      }
    }

    table[(diam+1)%2].clear();
  }

  booltable btable(GROUPSIZE);
  unsigned char *buf0 = new unsigned char[BUFSIZE];
  char *compressed_buf0 = new char[LZ4_compressBound(BUFSIZE)];
  unsigned char *buf1 = new unsigned char[BUFSIZE];
  char *compressed_buf1 = new char[LZ4_compressBound(BUFSIZE)];

  if (!resume) {
    for (hashset::iterator g = table[diam%2].begin(); g != table[diam%2].end(); g++)
      btable.set1(*g);
    table[diam%2].clear();
  }

  bool first = true;

  for (; count != GROUPSIZE; diam++) {
    char s0[100], s1[100];
    //    if (*s) unlink(s); // save disk space

    sprintf(s0, prefix, diam);
    sprintf(s1, prefix, diam-1);

    if (!resume)
      btable.save(s0); // make a copy of the data, to use it to fill in
    else {
      btable.load(s0); // initial run, load the data
      resume = false;
    }

    int f0 = open(s0, O_RDONLY), f1;
    if (f0 < 0) {
      std::cerr << "Error opening file " << s0 << " for reading\n";
      exit(1);
    }
    if (!first) {
      f1 = open(s1, O_RDONLY);
      if (f1 < 0) {
        std::cerr << "Error opening file " << s1 << " for reading\n";
        exit(1);
      }
    }

    for (activitymask g = 0; g < GROUPSIZE; g += 8*BUFSIZE) {
      int len;
      if (read (f0, &len, 4) != 4) {
	std::cerr << "Error reading from file " << s0 << " at position " << g << "\n";
	exit(1);
      }
      if (read (f0, compressed_buf0, len) != len) {
	std::cerr << "Error reading from file " << s0 << " at position " << g << "\n";
	exit(1);
      }
      LZ4_uncompress(compressed_buf0, (char *) buf0, BUFSIZE);

      if (first)
	memset (buf1, 0, BUFSIZE);
      else {
        if (read (f1, &len, 4) != 4) {
          std::cerr << "Error reading from file " << s1 << " at position " << g << "\n";
          exit(1);
        }
        if (read (f1, compressed_buf1, len) != len) {
          std::cerr << "Error reading from file " << s1 << " at position " << g << "\n";
          exit(1);
        }
        LZ4_uncompress(compressed_buf1, (char *) buf1, BUFSIZE);
      }
#pragma omp parallel for
      for (signed long gshift = 0; gshift < 8*BUFSIZE; gshift++) {
	if ((buf0[gshift/8] & pow2[gshift%8]) && !(buf1[gshift/8] & pow2[gshift%8])) {
	  treeaut t;
	  t.unpack(g+gshift);
	  btable.set1(t.pack());
	  btable.set1((t*A).pack());
	  btable.set1((t*B).pack());
	  btable.set1((t*C).pack());
	  btable.set1((t*D).pack());
	}
      }
    }      
    close(f0);
    if (first)
      first = false;
    else
      close(f1);

    oldcount = count;
    count = btable.size();
    std::cout << diam+1 << "\t" << count-oldcount << std::endl;
    output << diam+1 << "\t" << count-oldcount << std::endl;
    output.flush();
  }

  delete buf0;
  delete buf1;
  delete compressed_buf0;
  delete compressed_buf1;

  return 0;
}

#else // FULLTABLE

typedef std::set<htype> booltable;

#define MAXLEN 15

booltable table[MAXLEN][MAXLEN][MAXLEN][2];

int main(void) {
  table[0][0][0][0].insert(trivial.pack());

  for (int sum = 0; sum < 3*MAXLEN; sum++)
    for (int i = 0; i < sum && i < MAXLEN; i++)
      for (int j = 0; j < sum-i && j < MAXLEN; j++)
	for (int k = 0; k < sum-i-j && k < MAXLEN; k++)
	  for (int l = 0; l < 2; l++) {
	    booltable *thistable = table[i][j][k]+l, *lasttable, *nexttable;
	    if (i < MAXLEN-1) { // can multiply by a
	      if (i == 0)
		lasttable = NULL;
	      else
		lasttable = table[i-1][j][k]+(l+j+k)%2;
	      nexttable = table[i+1][j][k]+(l+j+k)%2;
	      for (booltable::iterator g = thistable->begin(); g != thistable->end(); g++) {
		treeaut t(*g);
		htype h = (t*A).pack();
		
		if (lasttable == NULL || lasttable->find(h) == lasttable->end())
		  nexttable->insert(h);
	      }
	    }
	    if (j < MAXLEN-1) { // can multiply by c
	      if (j == 0)
		lasttable = NULL;
	      else
		lasttable = table[i][j-1][k]+l;
	      nexttable = table[i][j+1][k]+l;
	      for (booltable::iterator g = thistable->begin(); g != thistable->end(); g++) {
		treeaut t(*g);
		htype h = (t*C).pack();
		
		if (lasttable == NULL || lasttable->find(h) == lasttable->end())
		  nexttable->insert(h);
	      }
	    }
	    if (k < MAXLEN-1) { // can multiply by d
	      if (k == 0)
		lasttable = NULL;
	      else
		lasttable = table[i][j][k-1]+l;
	      nexttable = table[i][j][k+1]+l;
	      for (booltable::iterator g = thistable->begin(); g != thistable->end(); g++) {
		treeaut t(*g);
		htype h = (t*D).pack();
		
		if (lasttable == NULL || lasttable->find(h) == lasttable->end())
		  nexttable->insert(h);
	      }
	    }
	  }

  for (int i = 0; i < MAXLEN; i++)
    for (int j = 0; j < MAXLEN; j++)
      for (int k = 0; k < MAXLEN; k++)
	for (int l = 0; l < 2; l++) {
	  long c = table[i][j][k][l].size();
	  if (c != 0)
	    std::cout << "table[" << i << "," << j << "," << k << "," << l << "] = " << c << std::endl;
	}
}
#endif
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
 *	216,356,488,772,1054,1660,2240,3448,4642,7128,
 *	9518,14392,19186,28984,38237,56452,74422,109676,143511,209636,
 *	274674,400028,518950,744680,966777,1385720,1787250,2540904,3278736,4642488,
 *	5931336,8292968,10596474,14788072,18777713,25981504,32988476,45487952,57193926,77818672,
 *	97700704,132478656,165252228,222000492,276453124,369614788,456758190,604793332,746199870,984109132,
 *	1206555972,1578342020,1931851085,2515791468,3038141194,3879710784,4672299320,5942511712,7086947781,8884186220,
 *	10540288059,13110253796,15378627270,18830675056,21957408389,26646407040,30723624604,36699612900,41999722934,49645170460,
 *	56023629761,64917048572,72613741621,83128358756,91647795632,102859681448,112207122317,124144579256,133181603189,143975261456,
 *	152527161217,162235208992,168832564233,175308584548,179754447175,183058600508,183726801446,181891660420,179514410319,173944013164,
 *	167600290329,157287440748,148460792785,135816646804,124622861278,109848342220,98381675585,84128387492,72767273549,59468201740,
 *	49982613608,39472047396,31768517757,23651391364,18338413444,13096559100,9594098205,6343409328,4435828834,2786268896,
 *	1789742612,1000052196,610961926,330522044,179247531,79889344,41344984,19518176,8461976,2847708,
 *	1257301,529652,144691,24152,12156,6456,1382]
 */
/* 1 4 8 24 56 136 [344] */

/* grigorchuk portrait:
p := function(g,n)
    local s, i, j;
    s := "[";
    for i in [n,n-1..1] do
        for j in [0,2^i..2^n-2^i] do
            if IsOddInt(QuoInt((j+1)^g-(j+1)+2^n,2^(i-1))) then Append(s,"@"); else Append(s,"."); fi;
            Append(s," ");
        od;
        Remove(s);
        Append(s,"|");
    od;
    Remove(s);
    Append(s,"]");
    return s;
end;
*/
