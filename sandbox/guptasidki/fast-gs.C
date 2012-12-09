#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <mt.h>

#define DEGREE 3
#define GPSIZE 3

#define LEVEL 4
#define SIZE (1+DEGREE+DEGREE*DEGREE+DEGREE*DEGREE*DEGREE)
#define LOGORDER 19
int keys[] = { 0,
	       1,                     14,
	       2,      6,             15,                     28,
	       3, 4,   7,     11,     16,17,  20,     24,     29,30, 33, 37 };

#define PERFECT_HASHING
#define LONGLONGHASH

#ifdef PERFECT_HASHING
typedef unsigned long htype;
#define HSIZE 1162261467
#elif defined (LONGLONGHASH)
typedef unsigned long long htype;
#define HSIZE 0x200000
#define HFREE 0xffffffffffffffffLL
#else
typedef unsigned long htype;
#define HSIZE 0x400000
#define HFREE 0xffffffff
#endif

int child[SIZE][DEGREE];
int rec_init (int i, int level)
{
  int pos = i+1;
  for (int d = 0; d < DEGREE; d++) {
    if (level == 1)
      child[i][d] = 0;
    else {
      child[i][d] = pos;
      pos = rec_init (pos, level-1);
    }
  }
  return pos;
}
const int product[GPSIZE][GPSIZE] = {
  { 0, 1, 2 },
  { 1, 2, 0 },
  { 2, 0, 1 } };
const int inverse[GPSIZE] =
  { 0, 2, 1 };
const int action[GPSIZE][DEGREE] = {
  { 0, 1, 2 },
  { 1, 2, 0 },
  { 2, 0, 1 } };

struct treeaut {
  char portrait[SIZE];

  treeaut () {}
  treeaut (char t) {
    if (!child[0][0]) rec_init (0, LEVEL);

    switch (t) {
    case '1':
      for (int i = 0; i < SIZE; i++) portrait[i] = 0;
      break;
    case 'a': case 'A':
      portrait[0] = 1 + (t == 'A');
      for (int i = 1; i < SIZE; i++) portrait[i] = 0;
      break;
    case 't': case 'T':
      for (int i = 0; i < SIZE; i++) portrait[i] = 0;
      for (int i = 0; child[i][0]; i = child[i][0]) {
	portrait[child[i][1]] = 1 + (t == 'T');
	portrait[child[i][2]] = 2 - (t == 'T');
      }
    }
  }
  int activity (void) { return portrait[0] == 0; }
  bool trivial (void) {
    for (int i = 0; i < SIZE; i++) if (portrait[i]) return false;
    return true;
  }
  htype perfect_packed (void) {
    htype w = 0;
    for (int i = 0; i < LOGORDER; i++)
	w = GPSIZE*w + portrait[keys[i]];
    return w;
  }
  htype packed (void) {
    htype w = 0;
    for (int i = 0; i < SIZE; i++)

	w = GPSIZE*w + portrait[i];
    return w;
  }
  void friend multiply_rec (char p[], int i, char q[], int j, char r[], int k)
    // multiply the subtree at q[j] with that at r[k], put result in p[i]
  {
    if (child[i][0]) {
      for (int d = 0; d < DEGREE; d++) {
	multiply_rec (p, child[i][d], q, child[j][d], r, child[k][action[q[j]][d]]);
      }
    }
    p[i] = product[q[j]][r[k]];
  }
  treeaut &operator *= (treeaut &x) {
    multiply_rec (portrait, 0, portrait, 0, x.portrait, 0);
    return *this;
  }
  treeaut operator * (treeaut &x) {
    treeaut n;

    multiply_rec (n.portrait, 0, portrait, 0, x.portrait, 0);
    return n;
  }
  void friend invert_rec (char p[], int i, char q[], int j)
    // invert the subtree at q[j], put result in p[i]
  {
    if (child[i][0]) {
      for (int d = 0; d < DEGREE; d++) {
	invert_rec (p, child[i][d], q, child[j][action[inverse[q[j]]][d]]);
      }
    }
    p[i] = inverse[q[j]];
  }
  treeaut operator -(void) {
    treeaut n;

    invert_rec (n.portrait, 0, portrait, 0);
    return n;
  }
};

void rec_print (std::ostream &s, char p[], int i)
  // print subtree at p[i]
{
  if (child[i][0]) {
    s << "[" << char('0'+p[i]);
    for (int d = 0; d < DEGREE; d++)
      s << " ", rec_print (s, p, child[i][d]);
    s << "]";
  } else
    s << char('0'+p[i]);
}
std::ostream &operator <<(std::ostream &s, treeaut t)
{
  rec_print (s, t.portrait, 0);
  return s;
}

#ifdef PERFECT_HASHING
char *hdata;
#else
htype hkey[HSIZE];
char hdata[HSIZE];
#endif

htype hashkey (treeaut &t)
{
#ifdef PERFECT_HASHING
  return t.perfect_packed();
#else
  //  htype p = t.packed();
  htype p = t.perfect_packed();

  unsigned long x = p & (HSIZE-1);

  while (hkey[x] != p) {
    if (hkey[x] == HFREE) { hkey[x] = p; return x; }
    x = (x + 1) & (HSIZE-1);
  }
  return x;
#endif
}

void print_key (htype t)
{
  char s[100];
  strcpy (s, "????????????????????????????????????????");
#ifdef PERFECT_HASHING
  for (int i = LOGORDER-1; i >= 0; i--) {
    s[keys[i]] = '0'+(t % GPSIZE);
    t /= GPSIZE;
  }
#else
  for (int i = SIZE-1; i >= 0; i--) {
    s[i] = '0'+(t % GPSIZE);
    t /= GPSIZE;
  }
#endif
  std::cout << s << "\n";
}

treeaut trivial('1'), gen_a('a'), gen_A('A'), gen_t('t'), gen_T('T');

int level, maxdepth, diameter;

void fillin (treeaut t, int depth)
{
  // for (int i = 0; i < depth; i++) std::cout << "  "; std::cout << depth << t << "\n";
  htype k = hashkey(t);
  if (hdata[k] == -1 || hdata[k] > depth) {
    hdata[k] = depth;
    // std::cout << t << " depth: " << depth << "\n";
    if (depth == maxdepth) return;
    fillin (t * gen_a, depth+1);
    fillin (t * gen_A, depth+1);
    fillin (t * gen_t, depth+1);
    fillin (t * gen_T, depth+1);
  }
}

void printgeod (treeaut t, int depth)
{
  static char mem[256];

  htype k = hashkey(t);
  if (hdata[k] != depth) return;

  if (depth == maxdepth) {
    mem[depth] = 0;
    std::cout << mem+level << " " << t << "\n";
    return;
  }
  mem[depth] = 'a'; printgeod (t * gen_a, depth+1);
  mem[depth] = 'A'; printgeod (t * gen_A, depth+1);
  mem[depth] = 't'; printgeod (t * gen_t, depth+1);
  mem[depth] = 'T'; printgeod (t * gen_T, depth+1);
}

int main (int argc, char *argv[]) {
#ifdef PERFECT_HASHING
  hdata = (char *) malloc (HSIZE);
  memset (hdata, -1, HSIZE);
#else
  memset (hdata, -1, sizeof hdata);
  memset (hkey, -1, sizeof hkey);
#endif

  if (argc != 2) {
    std::cerr << "Use: " << argv[0] << " <maxdepth>\n";
    return -1;
  }
  maxdepth = atoi(argv[1]);

  int gpsize = 0;

  for (level = maxdepth-10; level >= 0; level --) {
    fillin (trivial, level);
    
    int hist[256], total = 0;
    memset (hist, 0, sizeof hist);

    for (int i = 0; i < HSIZE; i++) {
      if (hdata[i] >= 0)
	hist[hdata[i]]++, total++;
    }
    for (int i = level; i < 256 && hist[i]; i++)
      std::cout << "Len " << i-level << ": " << hist[i] << "\n";
    std::cout << "Total: " << total << "\n";
    std::cout.flush();
    if (total == HSIZE)
      break;
  }

  int hist[256];
  memset (hist, 0, sizeof hist);

  for (int i = 0; i < HSIZE; i++) {
    if (hdata[i] >= 0)
      hist[hdata[i]]++;
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

//LEVEL 1: [1,2]
//LEVEL 2: [1,4,8,12,2]
//LEVEL 3: [1,4,8,16,32,
//          64,120,200,348,594,
//          604,168,28]
//LEVEL 4: [1,4,8,16,32,64,128,256,512,1020,
//          2020,4016,7920,15360,29760,
//          57136,108542,205856,389050,727212,
//          1346408,2466606,4465426,7948702,13882228,
//          23542792,38540412,60397600,89744530,124340300,
//          158279646,182281004,183685294,148746542,84882850,
//          29776588,5766398,578474,35610,4104,
//          1008,32
