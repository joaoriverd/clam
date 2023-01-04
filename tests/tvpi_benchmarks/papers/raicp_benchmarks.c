extern void __CRAB_assert(int);
extern void __CRAB_assume(int);
extern void __CRAB_get_range(double);
extern void __SEAHORN_error(int);

int rand();

/* The expected result is [0, 10] */
double bigLoop(double x, int N) {
  __CRAB_assume(x >= 0.0);
  __CRAB_assume(x <= 10.0);
  __CRAB_assume(N >= 1);
  __CRAB_assume(N <= 1000000);

  double a = 0.1;
  int i = 1;
  double y = x*x-x;
  __CRAB_get_range(y);

  if (y < 0.0) {
    __CRAB_get_range(x);
//    if (x > 1.2) { // tvpi is not accurate enough (yet) to prove that this cannot be the case.
    if (x > 4.0) {
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

/* The widening seems to be over-approximating this one. When increasing the widening delay to ~15 iterations,
 * the result seems to almost converge to [1.64234, 3.04445] before applying the widening. */
double sqrt1(double x) {
  __CRAB_assume(x >= 4.5);
  __CRAB_assume(x <= 5.5);

  double xn, xn1;

  xn = x/2.0;
  xn1 = 0.5*(xn + x/xn);
  __CRAB_get_range(xn1);

  while(xn-xn1 > 0.01) {
    xn = xn1;
    xn1 = 0.5*(xn + x/xn);
    __CRAB_get_range(xn1);
  }

  __CRAB_get_range(xn1);
  return xn1;
}

/* The widening seems to be over-approximating this one. When increasing the widening delay to ~15 iterations,
 * the result seems to almost converge to [1.00242, 7.48096] before applying the widening. */
double sqrt2(double x) {
  __CRAB_assume(x >= 5);
  __CRAB_assume(x <= 10);

  double xn, xn1;

  xn = x/2.0;
  xn1 = 0.5*(xn + x/xn);
  __CRAB_get_range(xn1);

  while(xn-xn1 > 0.01) {
    xn = xn1;
    xn1 = 0.5*(xn + x/xn);
  __CRAB_get_range(xn1);
  }

  __CRAB_get_range(xn1);
  return xn1;
}