extern void __CRAB_assert(int);
extern void __CRAB_assume(int);
extern void __CRAB_get_range(double);
extern void __SEAHORN_error(int);
extern int rand();

double artificial(double a, double b) {
  __CRAB_assume(a >= 0.0);
  __CRAB_assume(a <= 2.5);
  __CRAB_assume(b >= 0.0);
  __CRAB_assume(b <= 2.5);

  /* First part */
  double x;
  double sub = b - a;
  __CRAB_get_range(sub);
  if (sub >= 0.0) {
    x = b / (sub + 1.0);
  }
  else {
    x = 2*b;
  }
  __CRAB_get_range(x);

  /* Second part */
  double c = 0.5;
  double y = x*x-x;
  if (y < 0.0) {
    __CRAB_get_range(x);
    if (x > 4.0) {
      c = -2;
    }
  }

  __CRAB_get_range(c);
  while (rand()) {
    x = c * x;
  }

  __CRAB_get_range(x);
  return x;
}