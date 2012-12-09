#include <iostream>
#include <stdlib.h>

#define MAXLEN 1000
#define VALENCY 2

struct word {
  int len;
  char data[MAXLEN];

  word (void) { len = 0; }
  word (char *s) { len = strlen(s); strcpy (data, s); }
  word &operator *= (int i) {
    if (len > 0 && data[len-1] == i)
      len--;
    else if (len > 0 && data[len-1] > 'a' && i > 'a')
      data[len-1] = 'b'+'c'+'d'-data[len-1]-i;
    else data[len++] = i;
    return *this;
  }
  void reduce (void) {
    for (int i = 1; i < len; i++)
      if (data[i-1] == data[i]) {
	for (; i < len-1; i++) data[i-1] = data[i+1];
	len -= 2;
	reduce();
	return;
      }
  }
  char &operator [] (int i) { return data[i]; }
  int decompose (word w[VALENCY]) {
    for (int i = 0; i < VALENCY; i++)
      w[i].len = 0;
    int c = 0;
    for (int i = 0; i < len; i++)
      switch (data[i]) {
      case 'a':	c ^= 1;	break;
      case 'b': w[c] *= 'a'; w[c ^ 1] *= 'c'; break;
      case 'c': w[c] *= 'a'; w[c ^ 1] *= 'd'; break;
      case 'd': w[c ^ 1] *= 'b'; break;
    }
    return c;
  }
  bool check (bool fast = false) {
    if (0) {
      std::cout << "(";
      for (int i = 0; i < len; i++) std::cout << data[i];
      std::cout << ")";
    }
    if (len <= 3) return true;
    int state = 0x1000;
    for (int i = 0; i < len - fast; i++)
      switch (data[i] + state) {
      case 'a'+0x1000: state = 0x000; break;
      case 'b'+0x1000: state = 0x100; break;
      case 'c'+0x1000: state = 0x500; break;
      case 'd'+0x1000: state = 0x900; break;
					  
      case 'b'+0x000: state = 0x100; break;
      case 'a'+0x100: state = 0x200; break;
      case 'b'+0x200: state = 0x100; break;
      case 'c'+0x200: state = 0x500; break;
      case 'd'+0x200: state = 0x300; break;
      case 'a'+0x300: state = 0x400; break;
      case 'b'+0x400: return false;
      case 'c'+0x400: state = 0x500; break;
      case 'd'+0x400: return false;

      case 'c'+0x000: state = 0x500; break;
      case 'a'+0x500: state = 0x600; break;
      case 'b'+0x600: state = 0x100; break;
      case 'c'+0x600: state = 0x500; break;
      case 'd'+0x600: state = 0x700; break;
      case 'a'+0x700: state = 0x800; break;
      case 'b'+0x800: return false;
      case 'c'+0x800: return false;
      case 'd'+0x800: return false;

      case 'd'+0x000: state = 0x900; break;
      case 'a'+0x900: state = 0xa00; break;
      case 'd'+0xa00: return false;
      case 'b'+0xa00: state = 0x100; break;
      case 'c'+0xa00: state = 0x500; break;
      default:
        std::cerr << "INVALID STATE " << data[i] << " " << state << "\n";
      }
    if (1 || fast) {
      word w[VALENCY];
      decompose (w);
      for (int i = 0; i < VALENCY; i++)
        if (!w[i].check(fast)) return false;
    }
    return true;
  }
};

std::ostream &operator <<(std::ostream &s, word &w)
{
  for (int i = 0; i < w.len; i++)
    s << w.data[i];
  return s;
}

word w;
void search (int depth) {
  if (!w.check(true)) return;
  if (depth == 0) {
    if (w.check(false)) std::cout << w << "\n";
    // std::cout << w << ": " << w.check() << "\n";
    return;
  }
  w[w.len++] = 'a';
  w[w.len++] = 'b'; search (depth-1);
  w[w.len-1] = 'c'; search (depth-1);
  w[w.len-1] = 'd'; search (depth-1);
  w.len -= 2;
}  

int main (int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Use: " << argv[0] << " <word>\n";
    return -1;
  }
  //  word w(argv[1]);
  //  std::cout << "`" << w << "': " << w.check() << "\n";

  search (atoi(argv[1]));
}
