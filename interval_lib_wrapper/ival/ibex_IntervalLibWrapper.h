#ifndef _IBEX_INTERVALLIBWRAPPER_H_
#define _IBEX_INTERVALLIBWRAPPER_H_

#include <math.h>
#include <ival.h>
#include <ival-scalar.h>

#define IBEX_INTERVAL_LIB_NEG_INFINITY (-INFINITY)
#define IBEX_INTERVAL_LIB_POS_INFINITY (INFINITY)

/* An interval as ival writes it: its two ends, [NaN, NaN] being empty, as Gaol keeps it (16 bytes; a flag beside
   them made it 24, and IBEX's interval matrices are large). */
class IVAL_INTERVAL
{
	public:
		double inf;
		double sup;

	IVAL_INTERVAL(void) : inf(NAN), sup(NAN) {}
	IVAL_INTERVAL(double a, double b)
	{
		if (!(a <= b) || a == INFINITY || b == -INFINITY) inf = sup = NAN;
		else { inf = a == 0 ? 0.0 : a; sup = b == 0 ? 0.0 : b; }
	}
	IVAL_INTERVAL(double a)
	{
		if (a != a || a == INFINITY || a == -INFINITY) inf = sup = NAN;
		else inf = sup = a == 0 ? 0.0 : a;
	}
	/* a result of ival's, taken as it is: ordered with +0 ends, or [NaN, NaN] */
	static IVAL_INTERVAL of(double l, double h) { IVAL_INTERVAL r; r.inf = l; r.sup = h; return r; }
};

namespace ibex {
  typedef IVAL_INTERVAL interval_type_wrapper;

  /* the Hausdorff distance, rounded up */
  static inline double _interval_distance_wrapper(const interval_type_wrapper &x1, const interval_type_wrapper &x2)
  {
    double a = x1.inf, b = x2.inf, c = x1.sup, d = x2.sup, l1, h1, l2, h2;
    ival1_sub(a, a, b, b, &l1, &h1);
    ival1_sub(c, c, d, d, &l2, &h2);
    double u = fabs(l1) > fabs(h1) ? fabs(l1) : fabs(h1), v = fabs(l2) > fabs(h2) ? fabs(l2) : fabs(h2);
    if (u != u || v != v) return INFINITY;   /* an infinite end on both: the ends differ by inf - inf */
    return u > v ? u : v;
  }
}

#endif /* _IBEX_INTERVALLIBWRAPPER_H_ */
