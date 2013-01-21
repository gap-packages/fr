#include <iostream.h>
#include <stdlib.h>
#include <mt.h>

// we work with a group G = <A,B>, where A acts on the top of the tree
// X^*, and B acts recursively as (partially) elements of A,B and
// (partially) permutations on the top.

// the model is: X=(Z/5)^2, A=Z/5 acts on the first coordinate, and
// b in B=Z/5 acts as [[b,b^-1,a,a^-1,1], b, b^-1, 1, 1]

#define MAXLEN 1000
#define VALENCY 5
#define GENORDER 5

#define TYPE_A 0x00
#define TYPE_B 0x10
#define TYPE_X 0x20
#define TYPE_Y 0x30

struct generator {
  int data;

  generator () {}
  generator (int t) {
    data = t + 1 + (ul_rand() % (GENORDER-1));
  }
  generator (int type, int value) { data = type + value; }
  int type() { return data & 0xf0; }
  int value() { return data & 0x0f; }
  bool trivial() { return value() == 0; }

  generator &operator += (int i) {
    data += i;
    if (value() >= GENORDER) data -= GENORDER;
  }
  generator &operator += (generator g) {
    data += g.value();
    if (value() >= GENORDER) data -= GENORDER;
  }
  generator operator + (int i) {
    int j = value() + i;
    if (j >= GENORDER) j -= GENORDER;
    return generator (type(), j);
  }
  generator operator -(void) {
    const static int INVERSE[] = { 0, 4, 3, 2, 1 };
    return generator(type(), INVERSE[value()]);
  }
  char name (void) {
    switch (type()) {
    case TYPE_A: return '0'+value();
    case TYPE_B: return 'a'+value();
    case TYPE_X: return 'A'+value();
    case TYPE_Y: return 'F'+value();
    }
  }
};

struct word {
  generator data[MAXLEN];
  int len;

  word (void) { len = 0; }
  word (int l) {
    len = l;
    for (int i = 0; i < len; i++)
      data[i] = generator(i & 1 ? TYPE_Y : TYPE_A);
  }
  word &operator *= (generator i) {
    if (len > 0 && data[len-1].type() == i.type()) {
      data[len-1] += i;
      if (data[len-1].trivial())
	len--;
    }
    else data[len++] = i;
    return *this;
  }
  void reduce (void) {
    for (int i = 1; i < len; i++)
      if (data[i-1].type() == data[i].type()) {
	data[i-1] += data[i];
	int skip = 1 + data[i-1].trivial();
	for (i++; i < len; i++) data[i-skip] = data[i];
	len -= skip;
	reduce();
	return;
      }
  }
  int decompose (word w[VALENCY], generator c) {
    for (int i = 0; i < len; i++)
      switch (data[i].type()) {
      case TYPE_A:
	c += data[i];
	break;
      case TYPE_B:
	w[c.value()] *= generator(TYPE_X,data[i].value());
	w[(c+1).value()] *= generator(TYPE_A,data[i].value());
	w[(c+2).value()] *= -generator(TYPE_A,data[i].value());
	break;
      case TYPE_X:
	w[c.value()] *= generator(TYPE_Y,data[i].value());
	w[(c+1).value()] *= generator(TYPE_A,data[i].value());
	w[(c+2).value()] *= -generator(TYPE_A,data[i].value());
	break;
      case TYPE_Y:
	w[c.value()] *= generator(TYPE_B,data[i].value());
	w[(c+1).value()] *= -generator(TYPE_B,data[i].value());
	w[(c+2).value()] *= generator(TYPE_A,data[i].value());
	w[(c+3).value()] *= -generator(TYPE_A,data[i].value());
	break;
      }
    return c.value();
  }
  word &subword(int p, int l) {
    word *w = new word;

    w->len = l;
    for (int i = 0; i < l; i++)
      w->data[i] = data[i+p];
    return *w;
  }
};

ostream &operator <<(ostream &s, word &w)
{
  for (int i = 0; i < w.len; i++)
    s << w.data[i].name();
  return s;
}

bool is_trivial (word &w)
{
  if (w.len == 0)
    return true;
  if (w.len == 1)
    return false;
  word d[VALENCY];
  int c = w.decompose(d,generator(TYPE_A,0));
  if (c != 0) return false;
  for (int i = 0; i < VALENCY; i++)
    if (!is_trivial(d[i]))
      return false;
  return true;
}

int five_order (word w, int level)
{
  if (w.len == 0)
    return 0;
  if (w.len == 1)
    return 1;

  if (level > 20) cerr << "!!!" << w << " " << level << "\n";

  word d[VALENCY];
  int c = w.decompose(d,generator(TYPE_A,0));

  if (c == 0) {
    int m = 0;
    for (int i = 0; i < VALENCY; i++) {
      int n = five_order(d[i],level+1);
      if (n > m) m = n;
    }
    return m;
  }
  for (int i = 1; i < VALENCY; i++)
    c = w.decompose(d,generator(TYPE_A,c));
  return 1+five_order(d[0],level+1);
}

int main (int argc, char *argv[]) {
  if (argc != 6) {
    cerr << "Use: " << argv[0] << " <totlen> <lena> <lenb> <num> <seed>\n";
    return -1;
  }
  s_rand (atoi(argv[5]));

  int totlen = atoi(argv[1]),
    lena = atoi(argv[2]),
    lenb = atoi(argv[3]),
    num = atoi(argv[4]);

#ifdef PRINT_DECOMPOSITION
  for (int i = 0; i < num; i++) {
    word w(totlen), d[VALENCY];
    int c;
    c = w.decompose(d,generator(TYPE_A,0));
    cout << "WORD `" << w << "' <";
    for (int j = 0; j < VALENCY; j++) {
      if (j) cout << ",";
      cout << d[j];
    }
    cout << ">" << c << "\n";
  }
#elif defined(PRINT_ORDER)
  for (int i = 0; i < num; i++) {
    word w(totlen);
    cout << "WORD `" << w << "': order 5^" << five_order(w,0) << "\n";
  }
#else
  int triv = 0;
  for (int i = 0; i < num; i++) {
    word w(totlen);
    if (is_trivial(w)) {
      triv++;
      cout << "TRIVIAL `" << w << "'\n";
    }
  }
  cout << "num: " << num << ", trivial: " << triv << "\n";
#endif
}
