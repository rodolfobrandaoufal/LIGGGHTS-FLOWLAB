/* ----------------------------------------------------------------------
    LIGGGHTS - DEM simulation engine
    LIGGGHTS is open-source, distributed under the terms of the GNU Public
    License, version 2 or later.

    Contributing author: LIGGGHTS modernization branch (tests/quick, P0-18)

    Unit test for Utils::AbstractFactory: the style-table key must hold the
    full 64-bit contact-model hash. Two variants whose hashes differ only in
    the upper 32 bits must be registered without a collision warning and
    looked up independently. With the legacy int key the second addStyle
    prints "Style collision" and both lookups return the second creator.
    Exit code 0 = pass.
------------------------------------------------------------------------- */

#include "utils.h"
#include <cstdio>
#include <sstream>

using namespace LIGGGHTS;

namespace {
  struct DummyParent {};
  struct IDummy {
    typedef DummyParent ParentType;
    int which;
    explicit IDummy(int w) : which(w) {}
  };

  IDummy * createA(LAMMPS_NS::LAMMPS *, DummyParent *, int64_t) { static IDummy a(1); return &a; }
  IDummy * createB(LAMMPS_NS::LAMMPS *, DummyParent *, int64_t) { static IDummy b(2); return &b; }

  class Factory : public Utils::AbstractFactory<IDummy> {
  public:
    Factory() {}
  };
}

int main()
{
  int fails = 0;

  // the key type must be able to carry every bit of the hash type
  Factory f;
  const int64_t lo = Utils::generate_gran_hashcode(1, 2, 3, 4, 5);
  const int64_t hi = lo + (static_cast<int64_t>(1) << 32);   // same low 32 bits
  if (static_cast<int>(lo) != static_cast<int>(hi)) { std::printf("test setup error\n"); return 2; }

  std::stringstream captured;
  std::streambuf *old = std::cerr.rdbuf(captured.rdbuf());
  f.addStyle("gran", lo, &createA);
  f.addStyle("gran", hi, &createB);
  std::cerr.rdbuf(old);

  if (captured.str().find("collision") != std::string::npos) {
    std::printf("FAIL: collision warning for hashes differing only in the upper 32 bits: %s", captured.str().c_str());
    fails++;
  }
  IDummy *a = f.create("gran", lo, NULL, NULL);
  IDummy *b = f.create("gran", hi, NULL, NULL);
  if (!a || a->which != 1) { std::printf("FAIL: lookup of the low hash\n"); fails++; }
  if (!b || b->which != 2) { std::printf("FAIL: lookup of the high hash\n"); fails++; }
  if (!f.hasStaticVariant("gran", hi) || !f.hasStaticVariant("gran", lo)) { std::printf("FAIL: hasStaticVariant\n"); fails++; }
  if (f.hasStaticVariant("gran", hi + 1)) { std::printf("FAIL: hasStaticVariant of an unregistered hash\n"); fails++; }

  std::printf(fails ? "P0-18 factory key: FAIL\n" : "P0-18 factory key: PASS\n");
  return fails ? 1 : 0;
}
