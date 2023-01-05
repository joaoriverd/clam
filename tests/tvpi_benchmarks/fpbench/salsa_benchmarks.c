#include <math.h>
#define TRUE 1
#define FALSE 0

extern void __CRAB_assert(int);
extern void __CRAB_assume(int);
extern void __CRAB_get_range(double);
extern void __SEAHORN_error(int);

double PID(double m, double kp, double ki, double kd, double c) {
  __CRAB_assume(m >= -4.2);
  __CRAB_assume(m <= -4.0);
  __CRAB_assume(c >=  0.1);
  __CRAB_assume(c <=  0.2);
  __CRAB_assume(kp >= 9.454);
  __CRAB_assume(kp <= 9.50);
  __CRAB_assume(ki >= 0.69006);
  __CRAB_assume(ki <= 0.8);
  __CRAB_assume(kd >= 2.8454);
  __CRAB_assume(kd <= 2.98);

  double dt = 0.5;
  double invdt = 2.0;
  double e;
  double p;
  double d;
  double r;
  double m_1 = m;
  double eold = 0.0;
  double i = 0.0;
  double t = 0.0;
  while (t < 100.0) {
    e = c - m_1;
    p = kp * e;
    i = i + ((ki * dt) * e);
    d = (kd * invdt) * (e - eold);
    r = (p + i) + d;
    m_1 = m_1 + (0.01 * r);
    __CRAB_get_range(m_1);
    eold = e;
    t = t + dt;
  }
  __CRAB_get_range(m_1);
  return m_1;
}

#if 0
double Odometry(double sr_42_, double sl_42_) {
  __CRAB_assume(sr_42_ >= 0.05);
  __CRAB_assume(sr_42_ <= 6.28318530718);
  __CRAB_assume(sl_42_ >= 0.05);
  __CRAB_assume(sl_42_ <= 6.28318530718);

	double inv_l = 0.1;
	double c = 12.34;
	double delta_dl = 0.0;
	double delta_dr = 0.0;
	double delta_d = 0.0;
	double delta_theta = 0.0;
	double arg = 0.0;
	double cosi = 0.0;
	double x = 0.0;
	double sini = 0.0;
	double y = 0.0;
	double theta = -0.985;
	double t = 0.0;
	double tmp = sl_42_;
	double sl = sl_42_;
	double sr = sr_42_;
	double j = 0.0;
	while (t < 1000.0) {
		delta_dl = c * sl;
		delta_dr = c * sr;
		delta_d = (delta_dl + delta_dr) * 0.5;
		delta_theta = (delta_dr - delta_dl) * inv_l;
		arg = theta + (delta_theta * 0.5);
		cosi = (1.0 - ((arg * arg) * 0.5)) + ((((arg * arg) * arg) * arg) * 0.0416666666);
		x = x + (delta_d * cosi);
    __CRAB_get_range(x);
		sini = (arg - (((arg * arg) * arg) * 0.1666666666)) + (((((arg * arg) * arg) * arg) * arg) * 0.008333333);
		y = y + (delta_d * sini);
		theta = theta + delta_theta;
		t = t + 1.0;
		tmp = sl;
		double tmp_2;
		if (j == 50.0) {
			tmp_2 = sr;
		} else {
			tmp_2 = sl;
		}
		sl = tmp_2;
		double tmp_3;
		if (j == 50.0) {
			tmp_3 = tmp;
		} else {
			tmp_3 = sr;
		}
		sr = tmp_3;
		double tmp_4;
		if (j == 50.0) {
			tmp_4 = 0.0;
		} else {
			tmp_4 = j + 1.0;
		}
		j = tmp_4;
	}
  __CRAB_get_range(x);
	return x;
}

double Runge_Kutta_4(double h, double y_n_42_, double c) {
	double sixieme = 1.0 / 6.0;
	double eps = 0.005;
	double k = 1.2;
	double y_n = y_n_42_;
	double i = 0.0;
	double e = 1.0;
	int tmp = e > eps;
	while (tmp) {
		double v = c - y_n;
		double k1 = (k * v) * v;
		double v_1 = c - (y_n + ((0.5 * h) * k1));
		double k2 = (k * v_1) * v_1;
		double v_2 = c - (y_n + ((0.5 * h) * k2));
		double k3 = (k * v_2) * v_2;
		double v_3 = c - (y_n + (h * k3));
		double k4 = (k * v_3) * v_3;
		double y_n_4 = y_n + ((sixieme * h) * (((k1 + (2.0 * k2)) + (2.0 * k3)) + k4));
		double i_5 = i + 1.0;
		double e_6 = e - eps;
		y_n = y_n_4;
		i = i_5;
		e = e_6;
		tmp = e > eps;
	}
	return fabs(e);
}

double Lead_lag_System(double y, double yd) {
	double eps = 0.01;
	double Dc = -1280.0;
	double Ac00 = 0.499;
	double Ac01 = -0.05;
	double Ac10 = 0.01;
	double Ac11 = 1.0;
	double Bc0 = 1.0;
	double Bc1 = 0.0;
	double Cc0 = 564.48;
	double Cc1 = 0.0;
	double yc = 0.0;
	double u = 0.0;
	double xc0 = 0.0;
	double xc1 = 0.0;
	double i = 0.0;
	double e = 1.0;
	int tmp = e > eps;
	while (tmp) {
		double v = y - yd;
		double tmp_1;
		if (v < -1.0) {
			tmp_1 = -1.0;
		} else if (1.0 < v) {
			tmp_1 = 1.0;
		} else {
			tmp_1 = v;
		}
		yc = tmp_1;
		u = (Cc0 * xc0) + ((Cc1 * xc1) + (Dc * yc));
		xc0 = (Ac00 * xc0) + ((Ac01 * xc1) + (Bc0 * yc));
		xc1 = (Ac10 * xc0) + ((Ac11 * xc1) + (Bc1 * yc));
		i = i + 1.0;
		e = fabs((yc - xc1));
		tmp = e > eps;
	}
	return xc1;
}

double Trapeze(double u) {
	double a = 0.25;
	double b = 5000.0;
	double n = 25.0;
	double h = (b - a) / n;
	double xb = 0.0;
	double r = 0.0;
	double xa = 0.25;
	int tmp = xa < 5000.0;
	while (tmp) {
		double v = xa + h;
		double tmp_1;
		if (v > 5000.0) {
			tmp_1 = 5000.0;
		} else {
			tmp_1 = v;
		}
		xb = tmp_1;
		double gxa = u / ((((((0.7 * xa) * xa) * xa) - ((0.6 * xa) * xa)) + (0.9 * xa)) - 0.2);
		double gxb = u / ((((((0.7 * xb) * xb) * xb) - ((0.6 * xb) * xb)) + (0.9 * xb)) - 0.2);
		r = r + (((gxa + gxb) * 0.5) * h);
		xa = xa + h;
		tmp = xa < 5000.0;
	}
	return r;
}

double Rocket_Trajectory(double Mf, double A) {
	double R = 6400000.0;
	double G = 6.67428e-11;
	double Mt = 5.9736e+24;
	double dt = 0.1;
	double T = 24.0 * 3600.0;
	double nombrepas = T / dt;
	double r0 = (400.0 * 10000.0) + R;
	double vr0 = 0.0;
	double teta0 = 0.0;
	double viss = sqrtf(((G * Mt) / r0));
	double vteta0 = viss / r0;
	double rf = R;
	double vrf = 0.0;
	double tetaf = 0.0;
	double vl = sqrtf(((G * Mt) / R));
	double vlrad = vl / r0;
	double vtetaf = 1.1 * vlrad;
	double t_i = 0.0;
	double mf_i = 0.0;
	double u1_i = 0.0;
	double u3_i = 0.0;
	double w1_i = 0.0;
	double w3_i = 0.0;
	double u2_i = 0.0;
	double u4_i = 0.0;
	double w2_i = 0.0;
	double w4_i = 0.0;
	double x = 0.0;
	double y = 0.0;
	double i = 1.0;
	double u1_im1 = r0;
	double u2_im1 = vr0;
	double u3_im1 = teta0;
	double u4_im1 = vteta0;
	double w1_im1 = rf;
	double w2_im1 = vrf;
	double w3_im1 = tetaf;
	double w4_im1 = vtetaf;
	double t_im1 = 0.0;
	double mf_im1 = Mf;
	int tmp = i < 2000000.0;
	while (tmp) {
		t_i = t_im1 + dt;
		mf_i = mf_im1 - (A * t_im1);
		u1_i = (u2_im1 * dt) + u1_im1;
		u3_i = (u4_im1 * dt) + u3_im1;
		w1_i = (w2_im1 * dt) + w1_im1;
		w3_i = (w4_im1 * dt) + w3_im1;
		u2_i = ((-G * (Mt / (u1_im1 * u1_im1))) * dt) + ((u1_im1 * u4_im1) * (u4_im1 * dt));
		u4_i = ((-2.0 * (u2_im1 * (u4_im1 / u1_im1))) * dt) + u4_im1;
		double tmp_1;
		if (mf_im1 > 0.0) {
			tmp_1 = ((A * w2_im1) / (Mf - (A * t_im1))) * dt;
		} else {
			tmp_1 = 0.0;
		}
		w2_i = (((-G * (Mt / (w1_im1 * w1_im1))) * dt) + ((w1_im1 * w4_im1) * (w4_im1 * dt))) + (tmp_1 + w2_im1);
		double tmp_2;
		if (mf_im1 > 0.0) {
			tmp_2 = A * ((w4_im1 / (Mf - (A * t_im1))) * dt);
		} else {
			tmp_2 = 0.0;
		}
		w4_i = ((-2.0 * (w2_im1 * (w4_im1 / w1_im1))) * dt) + (tmp_2 + w4_im1);
		x = u1_i * cosf(u3_i);
		y = u1_i * sinf(u3_i);
		i = i + 1.0;
		u1_im1 = u1_i;
		u2_im1 = u2_i;
		u3_im1 = u3_i;
		u4_im1 = u4_i;
		w1_im1 = w1_i;
		w2_im1 = w2_i;
		w3_im1 = w3_i;
		w4_im1 = w4_i;
		t_im1 = t_i;
		mf_im1 = mf_i;
		tmp = i < 2000000.0;
	}
	return x;
}

double Jacobi_Method(double a11, double a22, double a33, double a44, double b1, double b2, double b3, double b4) {
	double eps = 1e-17;
	double x_n1 = 0.0;
	double x_n2 = 0.0;
	double x_n3 = 0.0;
	double x_n4 = 0.0;
	double i = 0.0;
	double e = 1.0;
	double x1 = 0.0;
	double x2 = 0.0;
	double x3 = 0.0;
	double x4 = 0.0;
	int tmp = e > eps;
	while (tmp) {
		x_n1 = (((b1 / a11) - ((0.1 / a11) * x2)) - ((0.2 / a11) * x3)) + ((0.3 / a11) * x4);
		x_n2 = (((b2 / a22) - ((0.3 / a22) * x1)) + ((0.1 / a22) * x3)) - ((0.2 / a22) * x4);
		x_n3 = (((b3 / a33) - ((0.2 / a33) * x1)) + ((0.3 / a33) * x2)) - ((0.1 / a33) * x4);
		x_n4 = (((b4 / a44) + ((0.1 / a44) * x1)) - ((0.2 / a44) * x2)) - ((0.3 / a44) * x3);
		i = i + 1.0;
		e = fabsf((x_n4 - x4));
		x1 = x_n1;
		x2 = x_n2;
		x3 = x_n3;
		x4 = x_n4;
		tmp = e > eps;
	}
	return x2;
}

double Newton_Raphson_Method(double x0) {
	double eps = 0.0005;
	double x_n = 0.0;
	double e = 1.0;
	double x = 0.0;
	double i = 0.0;
	int tmp = (e > eps) && (i < 100000.0);
	while (tmp) {
		double f = ((((((x * x) * ((x * x) * x)) - ((10.0 * x) * ((x * x) * x))) + ((40.0 * x) * (x * x))) - ((80.0 * x) * x)) + (80.0 * x)) - 32.0;
		double ff = (((((5.0 * x) * ((x * x) * x)) - ((40.0 * x) * (x * x))) + ((120.0 * x) * x)) - (160.0 * x)) + 80.0;
		x_n = x - (f / ff);
		e = fabs((x - x_n));
		x = x_n;
		i = i + 1.0;
		tmp = (e > eps) && (i < 100000.0);
	}
	return x;
}

double Eigenvalue_Computation(double a11, double a12, double a13, double a14, double a21, double a22, double a23, double a24, double a31, double a32, double a33, double a34, double a41, double a42, double a43, double a44, double v1, double v2, double v3, double v4) {
	double eps = 0.0005;
	double vx = 0.0;
	double vy = 0.0;
	double vz = 0.0;
	double vw = 0.0;
	double i = 0.0;
	double v1_1 = v1;
	double v2_2 = v2;
	double v3_3 = v3;
	double v4_4 = v4;
	double e = 1.0;
	int tmp = e > eps;
	while (tmp) {
		vx = ((a11 * v1_1) + (a12 * v2_2)) + ((a13 * v3_3) + (a14 * v4_4));
		vy = ((a21 * v1_1) + (a22 * v2_2)) + ((a23 * v3_3) + (a24 * v4_4));
		vz = ((a31 * v1_1) + (a32 * v2_2)) + ((a33 * v3_3) + (a34 * v4_4));
		vw = ((a41 * v1_1) + (a42 * v2_2)) + ((a43 * v3_3) + (a44 * v4_4));
		i = i + 1.0;
		v1_1 = vx / vw;
		v2_2 = vy / vw;
		v3_3 = vz / vw;
		v4_4 = 1.0;
		e = fabs((1.0 - v1_1));
		tmp = e > eps;
	}
	return v1_1;
}

double Iterative_Gram_Schmidt(double Q11, double Q12, double Q13, double Q21, double Q22, double Q23, double Q31, double Q32, double Q33) {
	double eps = 5e-6;
	double h1 = 0.0;
	double h2 = 0.0;
	double h3 = 0.0;
	double qj1 = Q31;
	double qj2 = Q32;
	double qj3 = Q33;
	double r1 = 0.0;
	double r2 = 0.0;
	double r3 = 0.0;
	double r = ((qj1 * qj1) + (qj2 * qj2)) + (qj3 * qj3);
	double rjj = 0.0;
	double e = 10.0;
	double i = 1.0;
	double rold = sqrtf(r);
	int tmp = e > eps;
	while (tmp) {
		h1 = ((Q11 * qj1) + (Q21 * qj2)) + (Q31 * qj3);
		h2 = ((Q12 * qj1) + (Q22 * qj2)) + (Q32 * qj3);
		h3 = ((Q13 * qj1) + (Q23 * qj2)) + (Q33 * qj3);
		qj1 = qj1 - (((Q11 * h1) + (Q12 * h2)) + (Q13 * h3));
		qj2 = qj2 - (((Q21 * h1) + (Q22 * h2)) + (Q23 * h3));
		qj3 = qj3 - (((Q31 * h1) + (Q32 * h2)) + (Q33 * h3));
		r1 = r1 + h1;
		r2 = r2 + h2;
		r3 = r3 + h3;
		r = ((qj1 * qj1) + (qj2 * qj2)) + (qj3 * qj3);
		rjj = sqrtf(r);
		e = fabsf((1.0 - (rjj / rold)));
		i = i + 1.0;
		rold = rjj;
		tmp = e > eps;
	}
	return qj1;
}

#endif