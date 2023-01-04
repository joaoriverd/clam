extern void __CRAB_assert(int);
extern void __CRAB_assume(int);
extern void __CRAB_get_range(double);
extern void __SEAHORN_error(int);

double simple_conditional(double a, double b) {
  __CRAB_assume(a >= 0.0);
  __CRAB_assume(a <= 100.0);
  __CRAB_assume(b >= 0.0);
  __CRAB_assume(b <= 100.0);

  double r;

  if (b >= a) {
    r = b / (b - a + 0.5);
    __CRAB_get_range(r);
  }
  else {
    r = b / 0.5;
    __CRAB_get_range(r);
  }

  __CRAB_get_range(r);
  return r;
}

double simple_conditional_adjusted(double a, double b) {
  __CRAB_assume(a >= 0.0);
  __CRAB_assume(a <= 100.0);
  __CRAB_assume(b >= 0.0);
  __CRAB_assume(b <= 100.0);

  double r;

  double sub = b - a;
  if (sub >= 0.0) {
    __CRAB_get_range(sub);
    r = b / (sub + 0.5);
    __CRAB_get_range(r);
  }
  else {
    r = b / 0.5;
    __CRAB_get_range(r);
  }

  __CRAB_get_range(r);
  return r;
}