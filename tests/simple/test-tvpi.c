// RUN: %clam -O0 --crab-dom=oct --crab-check=assert --crab-sanity-checks "%s" 2>&1 | OutputCheck %s
// CHECK: ^1  Number of total safe checks$
// CHECK: ^0  Number of total error checks$
// CHECK: ^0  Number of total warning checks$

extern void __CRAB_assert(int);
extern void __CRAB_assume(int);
extern void __SEAHORN_error(int);

//int foo1(int x, int y) {
//  int z;
//  __CRAB_assume(x <= 10);
//  __CRAB_assume(x >= -10);
//  __CRAB_assume(y <= 5);
//  __CRAB_assume(y >= -5);
//
//  if (x > 0) {
//    z = x + y;
//  }
//  else {
//    z = x - y;
//  }
//
//  __CRAB_assert(z <= 100);
//  return z;
//}


double foo2(double x, double y) {
    double z1, z2;
    __CRAB_assume(x <= 12);
    __CRAB_assume(x >= 10);
    __CRAB_assume(y <= 5);
    __CRAB_assume(y >= 4);

    z1 = x / y;
    z2 = y / x;

    __CRAB_assert(z1 + z2 <= 100);
    return z1;
}
