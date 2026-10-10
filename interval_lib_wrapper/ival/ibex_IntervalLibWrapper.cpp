#include "ibex_Interval.h"

namespace ibex {

/* pi rounded down and up, and pi/2 and 2pi (exact scalings of them) */
static const double IVAL_PI_DN = 3.141592653589793116, IVAL_PI_UP = 3.141592653589793560;

/* Rounding upward from the start, as Gaol's init sets it: IBEX's own code takes upward as its default mode and returns
   to it ("fpu_round_up(); // default mode"). ival is tight in any mode; in this one ival-scalar.h bounds a sum or a
   product with one operation each, where to nearest it needs an error-free transformation. */
static int ival_upward __attribute__((unused)) = (fesetround(FE_UPWARD), 0);

// *** Deprecated ***
const Interval Interval::EMPTY_SET((IVAL_INTERVAL()));
const Interval Interval::ALL_REALS(-INFINITY, INFINITY);
const Interval Interval::ZERO(0.0);
const Interval Interval::ONE(1.0);
const Interval Interval::POS_REALS(0.0, INFINITY);
const Interval Interval::NEG_REALS(-INFINITY, 0.0);
const Interval Interval::PI(IVAL_PI_DN, IVAL_PI_UP);
const Interval Interval::TWO_PI(2.0 * IVAL_PI_DN, 2.0 * IVAL_PI_UP);
const Interval Interval::HALF_PI(0.5 * IVAL_PI_DN, 0.5 * IVAL_PI_UP);

const Interval& Interval::empty_set() {
	static Interval _empty_set((IVAL_INTERVAL()));
	return _empty_set;
}
const Interval& Interval::all_reals() {
	static Interval _all_reals(-INFINITY, INFINITY);
	return _all_reals;
}
const Interval& Interval::zero() {
	static Interval _zero(0.0);
	return _zero;
}
const Interval& Interval::one() {
	static Interval _one(1.0);
	return _one;
}
const Interval& Interval::pos_reals() {
	static Interval _pos_reals(0.0, INFINITY);
	return _pos_reals;
}
const Interval& Interval::neg_reals() {
	static Interval _neg_reals(-INFINITY, 0.0);
	return _neg_reals;
}
const Interval& Interval::pi() {
	static Interval _pi(IVAL_PI_DN, IVAL_PI_UP);
	return _pi;
}
const Interval& Interval::two_pi() {
	static Interval _two_pi(2.0 * IVAL_PI_DN, 2.0 * IVAL_PI_UP);
	return _two_pi;
}
const Interval& Interval::half_pi() {
	static Interval _half_pi(0.5 * IVAL_PI_DN, 0.5 * IVAL_PI_UP);
	return _half_pi;
}

std::ostream& operator<<(std::ostream& os, const Interval& x) {
	if (x.is_empty()) return os << "[ empty ]";
	return os << "[" << x.lb() << ", " << x.ub() << "]";
}

} // end namespace
