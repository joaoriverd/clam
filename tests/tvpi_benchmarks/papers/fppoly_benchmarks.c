// RUN: %clam -O0 --crab-dom=oct --crab-check=assert --crab-sanity-checks "%s" 2>&1 | OutputCheck %s
// CHECK: ^1  Number of total safe checks$
// CHECK: ^0  Number of total error checks$
// CHECK: ^0  Number of total warning checks$

extern void __CRAB_assert(int);
extern void __CRAB_assume(int);
extern void __CRAB_get_range(double);
extern void __SEAHORN_error(int);

#define M 150

double FPPoly_example(double X, double Y, double D) {
  /* Result:
   * TVPI: [-150, 150] // It can possibly be improved
   */
  __CRAB_assume(Y <=  M);
  __CRAB_assume(Y >= -M);
  __CRAB_assume(X <=  128);
  __CRAB_assume(X >= -128);
  __CRAB_assume(D <=  16);
  __CRAB_assume(D >=  1);

  double S, R;

  while (rand()) {
    S = Y;
    R = X - S;
    Y = X;

    if (R <= -D) {
      Y = S - D;
    }
    else if (D <= R) {
      Y = S + D;
    }
  }

  __CRAB_get_range(Y);
  return Y;
}