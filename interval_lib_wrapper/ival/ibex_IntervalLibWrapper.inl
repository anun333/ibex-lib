/* The ival backend (crnumerics' ival, 2026-10-09). The arithmetic is ival-scalar.h's, inline and tight, upward
   (IBEX's default mode, as Gaol's) or to nearest; every elementary function, pow, root and backward operator is
   ival's, tight (CORE-MATH, correctly rounded both ways), on x86-64 and aarch64 alike. Conventions follow the gaol
   wrapper's: an infinite double operand gives the empty set, log of an interval with no positive point is empty;
   the empty set is [NaN, NaN], as Gaol's. */
#ifndef _IBEX_INTERVALLIBWRAPPER_INL_
#define _IBEX_INTERVALLIBWRAPPER_INL_

#include <fenv.h>
#include "ibex_Exception.h"

namespace ibex {

inline void fpu_round_down() { fesetround(FE_DOWNWARD); }
inline void fpu_round_up() { fesetround(FE_UPWARD); }
inline void fpu_round_near() { fesetround(FE_TONEAREST); }

inline double previous_float(double x) { return nextafter(x, -INFINITY); }
inline double next_float(double x) { return nextafter(x, INFINITY); }

namespace ival_wrapper {
	/* [l, h] from ival: valid, or [NaN, NaN] for empty, as the type keeps it */
	inline IVAL_INTERVAL mk(double l, double h) { return IVAL_INTERVAL::of(l, h); }
	inline bool inf_d(double d) { return d == INFINITY || d == -INFINITY; }
	typedef void (*fn1)(const double *, const double *, double *, double *, size_t);
	inline IVAL_INTERVAL f1(fn1 f, const Interval& x) {
		if (x.is_empty()) return IVAL_INTERVAL();
		double a = x.lb(), b = x.ub(), l, h;
		f(&a, &b, &l, &h, 1);
		return mk(l, h);
	}
	/* The four operations: ival1's take [NaN, NaN] as empty and give it back, and a double operand d as the point
	   [d, d], which is empty for an infinite d or a NaN, as this wrapper's convention wants; so no check here. */
	inline IVAL_INTERVAL add(double al, double ah, double bl, double bh) { double l, h; ival1_add(al, ah, bl, bh, &l, &h); return mk(l, h); }
	inline IVAL_INTERVAL sub(double al, double ah, double bl, double bh) { double l, h; ival1_sub(al, ah, bl, bh, &l, &h); return mk(l, h); }
	inline IVAL_INTERVAL mul(double al, double ah, double bl, double bh) { double l, h; ival1_mul(al, ah, bl, bh, &l, &h); return mk(l, h); }
	inline IVAL_INTERVAL div(double al, double ah, double bl, double bh) { double l, h; ival1_div(al, ah, bl, bh, &l, &h); return mk(l, h); }
	inline IVAL_INTERVAL scale(double d, double bl, double bh) { double l, h; ival1_scale(d, bl, bh, &l, &h); return mk(l, h); }
	/* the hull of {x in X : f(x) in Y}, by an ival reverse operation taking (C, X) */
	typedef void (*fnr)(const double *, const double *, const double *, const double *, double *, double *, size_t);
	inline bool rev(fnr f, const Interval& y, Interval& x) {
		if (y.is_empty() || x.is_empty()) { x.set_empty(); return false; }
		double cl = y.lb(), ch = y.ub(), xl = x.lb(), xh = x.ub(), l, h;
		f(&cl, &ch, &xl, &xh, &l, &h, 1);
		x = Interval(mk(l, h));
		return !x.is_empty();
	}
}

inline Interval::Interval(const IVAL_INTERVAL& x) : itv(x) {
}

inline Interval& Interval::operator=(const IVAL_INTERVAL& x) {
	this->itv = x;
	return *this;
}

inline Interval& Interval::operator+=(double d) { itv = ival_wrapper::add(lb(), ub(), d, d); return *this; }
inline Interval& Interval::operator-=(double d) { itv = ival_wrapper::sub(lb(), ub(), d, d); return *this; }
inline Interval& Interval::operator*=(double d) { itv = ival_wrapper::scale(d, lb(), ub()); return *this; }
inline Interval& Interval::operator/=(double d) { itv = ival_wrapper::div(lb(), ub(), d, d); return *this; }
inline Interval& Interval::operator+=(const Interval& x) { itv = ival_wrapper::add(lb(), ub(), x.lb(), x.ub()); return *this; }
inline Interval& Interval::operator-=(const Interval& x) { itv = ival_wrapper::sub(lb(), ub(), x.lb(), x.ub()); return *this; }
inline Interval& Interval::operator*=(const Interval& x) { itv = ival_wrapper::mul(lb(), ub(), x.lb(), x.ub()); return *this; }
inline Interval& Interval::operator/=(const Interval& x) { itv = ival_wrapper::div(lb(), ub(), x.lb(), x.ub()); return *this; }

inline Interval Interval::operator-() const {
	if (is_empty()) return *this;
	return Interval(IVAL_INTERVAL(-ub(), -lb()));
}

inline Interval& Interval::div2_inter(const Interval& x, const Interval& y) {
	Interval out2;
	div2_inter(x, y, out2);
	return *this |= out2;
}

inline void Interval::set_empty() {
	*this = EMPTY_SET;
}

inline Interval& Interval::operator&=(const Interval& x) {
	if (is_empty()) return *this;
	if (x.is_empty() || lb() > x.ub() || x.lb() > ub()) { set_empty(); return *this; }
	itv = IVAL_INTERVAL(lb() > x.lb() ? lb() : x.lb(), ub() < x.ub() ? ub() : x.ub());
	return *this;
}

inline Interval& Interval::operator|=(const Interval& x) {
	if (is_empty()) { *this = x; return *this; }
	if (x.is_empty()) return *this;
	itv = IVAL_INTERVAL(lb() < x.lb() ? lb() : x.lb(), ub() > x.ub() ? ub() : x.ub());
	return *this;
}

inline double Interval::lb() const { return itv.inf; }
inline double Interval::ub() const { return itv.sup; }

inline double Interval::mid() const {   /* 1788's midpoint: to nearest; -DBL_MAX, 0 or DBL_MAX when unbounded */
	if (is_empty()) return NAN;
	double a = lb(), b = ub(), m;
	ival_mid(&a, &b, &m, 1);
	return m;
}

inline bool Interval::is_empty() const { return itv.inf != itv.inf; }
inline bool Interval::is_degenerated() const { return is_empty() || lb() == ub(); }
inline bool Interval::is_unbounded() const { return !is_empty() && (lb() == -INFINITY || ub() == INFINITY); }

inline double Interval::diam() const {   /* the width rounded up; 0 for the empty set, as IBEX's direct backend */
	if (is_empty()) return 0;
	double a = lb(), b = ub(), w;
	ival_wid(&a, &b, &w, 1);
	return w;
}

inline double Interval::mig() const {
	if (lb() > 0) return lb();
	else if (ub() < 0) return -ub();
	else return 0;
}

inline double Interval::mag() const {
	return fabs(lb()) > fabs(ub()) ? fabs(lb()) : fabs(ub());
}

inline Interval operator&(const Interval& x1, const Interval& x2) { Interval r(x1); r &= x2; return r; }
inline Interval operator|(const Interval& x1, const Interval& x2) { Interval r(x1); r |= x2; return r; }

inline Interval operator+(const Interval& x, double d) { Interval r(x); r += d; return r; }
inline Interval operator-(const Interval& x, double d) { Interval r(x); r -= d; return r; }
inline Interval operator*(const Interval& x, double d) { return Interval(ival_wrapper::scale(d, x.lb(), x.ub())); }
inline Interval operator/(const Interval& x, double d) { Interval r(x); r /= d; return r; }
inline Interval operator+(double d, const Interval& x) { return Interval(ival_wrapper::add(d, d, x.lb(), x.ub())); }
inline Interval operator-(double d, const Interval& x) { return Interval(ival_wrapper::sub(d, d, x.lb(), x.ub())); }
inline Interval operator*(double d, const Interval& x) { return Interval(ival_wrapper::scale(d, x.lb(), x.ub())); }
inline Interval operator/(double d, const Interval& x) { return Interval(ival_wrapper::div(d, d, x.lb(), x.ub())); }
inline Interval operator+(const Interval& x1, const Interval& x2) { return Interval(ival_wrapper::add(x1.lb(), x1.ub(), x2.lb(), x2.ub())); }
inline Interval operator-(const Interval& x1, const Interval& x2) { return Interval(ival_wrapper::sub(x1.lb(), x1.ub(), x2.lb(), x2.ub())); }
inline Interval operator*(const Interval& x1, const Interval& x2) { return Interval(ival_wrapper::mul(x1.lb(), x1.ub(), x2.lb(), x2.ub())); }
inline Interval operator/(const Interval& x1, const Interval& x2) { return Interval(ival_wrapper::div(x1.lb(), x1.ub(), x2.lb(), x2.ub())); }

inline Interval sqr(const Interval& x) {
	if (x.is_empty()) return x;
	double l, h; ival1_sqr(x.lb(), x.ub(), &l, &h);
	return Interval(ival_wrapper::mk(l, h));
}

inline Interval sqrt(const Interval& x) {
	if (x.is_empty()) return x;
	double l, h; ival1_sqrt(x.lb(), x.ub(), &l, &h);
	return Interval(ival_wrapper::mk(l, h));
}

inline Interval pow(const Interval& x, int n) {
	if (x.is_empty()) return x;
	double a = x.lb(), b = x.ub(), l, h;
	ival_pown(&a, &b, &n, &l, &h, 1);
	return Interval(ival_wrapper::mk(l, h));
}

inline Interval pow(const Interval& x, double d) {
	if (ival_wrapper::inf_d(d) || x.is_empty()) return Interval::empty_set();
	double a = x.lb(), b = x.ub(), l, h;
	ival_pow(&a, &b, &d, &d, &l, &h, 1);
	return Interval(ival_wrapper::mk(l, h));
}

inline Interval pow(const Interval& x, const Interval& y) {
	if (x.is_empty() || y.is_empty()) return Interval::empty_set();
	double a = x.lb(), b = x.ub(), c = y.lb(), d = y.ub(), l, h;
	ival_pow(&a, &b, &c, &d, &l, &h, 1);
	return Interval(ival_wrapper::mk(l, h));
}

/* the real n-th root: over the whole line for odd n, over x >= 0 for even n; a negative n is the reciprocal */
inline Interval root(const Interval& x, int n) {
	if (x.is_empty()) return x;
	double a = x.lb(), b = x.ub(), l, h;
	ival_rootn(&a, &b, &n, &l, &h, 1);
	return Interval(ival_wrapper::mk(l, h));
}

inline Interval exp(const Interval& x) { return Interval(ival_wrapper::f1(ival_exp, x)); }

inline Interval log(const Interval& x) {
	if (x.is_empty() || x.ub() <= 0) return Interval::empty_set();
	return Interval(ival_wrapper::f1(ival_log, x));
}

inline Interval cos(const Interval& x) { return Interval(ival_wrapper::f1(ival_cos, x)); }
inline Interval sin(const Interval& x) { return Interval(ival_wrapper::f1(ival_sin, x)); }
inline Interval tan(const Interval& x) { return Interval(ival_wrapper::f1(ival_tan, x)); }
inline Interval acos(const Interval& x) { return Interval(ival_wrapper::f1(ival_acos, x)); }
inline Interval asin(const Interval& x) { return Interval(ival_wrapper::f1(ival_asin, x)); }
inline Interval atan(const Interval& x) { return Interval(ival_wrapper::f1(ival_atan, x)); }
inline Interval cosh(const Interval& x) { return Interval(ival_wrapper::f1(ival_cosh, x)); }
inline Interval sinh(const Interval& x) { return Interval(ival_wrapper::f1(ival_sinh, x)); }
inline Interval tanh(const Interval& x) { return Interval(ival_wrapper::f1(ival_tanh, x)); }
inline Interval acosh(const Interval& x) { return Interval(ival_wrapper::f1(ival_acosh, x)); }
inline Interval asinh(const Interval& x) { return Interval(ival_wrapper::f1(ival_asinh, x)); }
inline Interval atanh(const Interval& x) { return Interval(ival_wrapper::f1(ival_atanh, x)); }

inline Interval abs(const Interval& x) { return Interval(ival_wrapper::f1(ival_abs, x)); }

inline Interval max(const Interval& x, const Interval& y) {
	if (x.is_empty() || y.is_empty()) return Interval::empty_set();
	return Interval(IVAL_INTERVAL(x.lb() > y.lb() ? x.lb() : y.lb(), x.ub() > y.ub() ? x.ub() : y.ub()));
}

inline Interval min(const Interval& x, const Interval& y) {
	if (x.is_empty() || y.is_empty()) return Interval::empty_set();
	return Interval(IVAL_INTERVAL(x.lb() < y.lb() ? x.lb() : y.lb(), x.ub() < y.ub() ? x.ub() : y.ub()));
}

inline Interval floor(const Interval& x) { return Interval(ival_wrapper::f1(ival_floor, x)); }
inline Interval ceil(const Interval& x) { return Interval(ival_wrapper::f1(ival_ceil, x)); }

/* the integers in x, as an interval: [ceil(lb), floor(ub)], empty when there are none */
inline Interval integer(const Interval& x) {
	if (x.is_empty()) return x;
	return Interval(IVAL_INTERVAL(::ceil(x.lb()), ::floor(x.ub())));
}

/* y = x1 * x2: x1 within {x1 : x1 * x2 in y for some x2 in X2}, then x2 likewise (1788's mulRev) */
inline bool bwd_mul(const Interval& y, Interval& x1, Interval& x2) {
	if (y.is_empty() || x1.is_empty() || x2.is_empty()) { x1.set_empty(); x2.set_empty(); return false; }
	double b0 = x2.lb(), b1 = x2.ub(), c0 = y.lb(), c1 = y.ub(), a0 = x1.lb(), a1 = x1.ub(), l, h;
	ival_mulrev(&b0, &b1, &c0, &c1, &a0, &a1, &l, &h, 1);
	x1 = Interval(ival_wrapper::mk(l, h));
	if (x1.is_empty()) { x2.set_empty(); return false; }
	b0 = x1.lb(); b1 = x1.ub(); a0 = x2.lb(); a1 = x2.ub();
	ival_mulrev(&b0, &b1, &c0, &c1, &a0, &a1, &l, &h, 1);
	x2 = Interval(ival_wrapper::mk(l, h));
	if (x2.is_empty()) { x1.set_empty(); return false; }
	return true;
}

inline bool bwd_sqr(const Interval& y, Interval& x) { return ival_wrapper::rev(ival_sqrrev, y, x); }

inline bool bwd_pow(const Interval& y, int expon, Interval& x) {
	if (y.is_empty() || x.is_empty()) { x.set_empty(); return false; }
	double cl = y.lb(), ch = y.ub(), xl = x.lb(), xh = x.ub(), l, h;
	ival_pownrev(&cl, &ch, &xl, &xh, &expon, &l, &h, 1);
	x = Interval(ival_wrapper::mk(l, h));
	return !x.is_empty();
}

/* y = x1^x2: x1 within powRev1, then x2 within powRev2 (1788) */
inline bool bwd_pow(const Interval& y, Interval& x1, Interval& x2) {
	if (y.is_empty() || x1.is_empty() || x2.is_empty()) { x1.set_empty(); x2.set_empty(); return false; }
	double c0 = y.lb(), c1 = y.ub(), b0 = x2.lb(), b1 = x2.ub(), a0 = x1.lb(), a1 = x1.ub(), l, h;
	ival_powrev1(&b0, &b1, &c0, &c1, &a0, &a1, &l, &h, 1);
	x1 = Interval(ival_wrapper::mk(l, h));
	if (x1.is_empty()) { x2.set_empty(); return false; }
	a0 = x1.lb(); a1 = x1.ub();
	ival_powrev2(&a0, &a1, &c0, &c1, &b0, &b1, &l, &h, 1);
	x2 = Interval(ival_wrapper::mk(l, h));
	if (x2.is_empty()) { x1.set_empty(); return false; }
	return true;
}

inline bool bwd_cos(const Interval& y, Interval& x) { return ival_wrapper::rev(ival_cosrev, y, x); }
inline bool bwd_sin(const Interval& y, Interval& x) { return ival_wrapper::rev(ival_sinrev, y, x); }
inline bool bwd_tan(const Interval& y, Interval& x) { return ival_wrapper::rev(ival_tanrev, y, x); }
inline bool bwd_cosh(const Interval& y, Interval& x) { return ival_wrapper::rev(ival_coshrev, y, x); }

inline bool bwd_sinh(const Interval& y, Interval& x) {   /* sinh is a bijection: x within asinh(y) */
	x &= asinh(y);
	return !x.is_empty();
}

inline bool bwd_tanh(const Interval& y, Interval& x) {   /* tanh is a bijection onto (-1, 1): x within atanh(y) */
	x &= atanh(y);
	return !x.is_empty();
}

inline bool bwd_abs(const Interval& y, Interval& x) { return ival_wrapper::rev(ival_absrev, y, x); }

} // end namespace ibex

#endif /* _IBEX_INTERVALLIBWRAPPER_INL_ */
