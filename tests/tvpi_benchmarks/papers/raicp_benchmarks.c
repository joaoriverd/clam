// RUN: %clam -O0 --crab-dom=oct --crab-check=assert --crab-sanity-checks "%s" 2>&1 | OutputCheck %s
// CHECK: ^1  Number of total safe checks$
// CHECK: ^0  Number of total error checks$
// CHECK: ^0  Number of total warning checks$

extern void __CRAB_assert(int);
extern void __CRAB_assume(int);
extern void __CRAB_get_range(double);
extern void __SEAHORN_error(int);

double bigLoop(double x, int N) {
  __CRAB_assume(x >= 0.0);
  __CRAB_assume(x <= 10.0);
  __CRAB_assume(N >= 1);
  __CRAB_assume(N <= 1000000);

  double a = 0.1;
  int i = 1;
  double y = x*x+x; // should be x*x-x?
  __CRAB_get_range(y);

  if (y < -0.1) {
    if (x > 1.2) {
      a = -2;
    }
  }

  while (N > i) {
    x = a * x;
    i = i + 1;
  }

  __CRAB_get_range(x);
  return x;
}

#if 0
double sqrt1(double x) {
  __CRAB_assume(x >= 4.5);
  __CRAB_assume(x <= 5.5);

  double xn, xn1;

  xn = x/2.0;
  xn1 = 0.5*(xn + x/xn);

  while(xn-xn1 > 0.01) {
    xn = xn1;
    xn1 = 0.5*(xn + x/xn);
  }

  return xn1;
}

double sqrt2(double x) {
  __CRAB_assume(x >= 5);
  __CRAB_assume(x <= 10);

  double xn, xn1;

  xn = x/2.0;
  xn1 = 0.5*(xn + x/xn);

  while(xn-xn1 > 0.01) {
    xn = xn1;
    xn1 = 0.5*(xn + x/xn);
  }

  return xn1;
}
#endif