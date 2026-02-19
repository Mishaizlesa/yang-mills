#include <cmath>
#include <cstdint>

using dd = struct DoubleDouble {
    double hi, lo;
};


inline void TwoSum(double a, double b, double& s, double& e) {
    s = a + b;
    double v = s - a;
    e = (a - (s - v)) + (b - v);
}

inline void TwoDiff(double a, double b, double& s, double& e) {
    s = a - b;
    double v = s - a;
    e = (a - (s - v)) - (b + v);
}

inline void TwoProduct(double a, double b, double& p, double& e) {
    p = a * b;
    constexpr double factor = (1LL << 27) + 1;  // 2²⁷ + 1 = 134217729
    double a1 = a * factor;  double a2 = a - a1;  a1 -= a2;  // a = a1 + a2
    double b1 = b * factor;  double b2 = b - b1;  b1 -= b2;  // b = b1 + b2
    e = a2*b2 - (((p - a1*b1) - a2*b1) - a1*b2);
}

inline DoubleDouble mul_small(const DoubleDouble& a, double b) {
    double ph, eh, pl, el;
    TwoProduct(a.hi, b, ph, eh);
    TwoProduct(a.lo, b, pl, el);
    double s, e;
    TwoSum(ph, pl, s, e);
    e += eh + el;
    DoubleDouble r{s, e};
    return r;
}

inline DoubleDouble add_dd(const DoubleDouble& a, const DoubleDouble& b) {
    double s1, s2, t1, t2;
    TwoSum(a.hi, b.hi, s1, s2);
    TwoSum(a.lo, b.lo, t1, t2);
    double s, e;
    TwoSum(s1, t1, s, e);
    e += s2 + t2;
    return {s, e};
}

inline DoubleDouble sub_dd(const DoubleDouble& a, const DoubleDouble& b) {
    return add_dd(a, {-b.hi, -b.lo});
}

inline DoubleDouble mul_dd(const DoubleDouble& a, const DoubleDouble& b) {
    double ph, eh1, eh2, pl, el;
    TwoProduct(a.hi, b.hi, ph, eh1);
    TwoProduct(a.hi, b.lo, pl, el);  eh2 = eh1 + pl;
    TwoProduct(a.lo, b.hi, pl, el);  eh2 += pl;
    TwoProduct(a.lo, b.lo, pl, el);
    double s, e;
    TwoSum(ph, eh2, s, e);
    e += el;
    return {s, e};
}

inline DoubleDouble to_dd(double x) {
    return {x, 0.0};
}

inline DoubleDouble operator-(const DoubleDouble& a) {
    return {-a.hi, -a.lo};
}

inline DoubleDouble reciprocal_dd(const DoubleDouble& a) {
    if (a.hi == 0.0) {
        return {0.0 / 0.0, 0.0 / 0.0}; 
    }

    double q = 1.0 / a.hi;

    double ph, eh;
    TwoProduct(a.hi, q, ph, eh);
    DoubleDouble aq = {ph, eh + a.lo * q};

    DoubleDouble one = {1.0, 0.0};
    DoubleDouble r = sub_dd(one, aq);

    DoubleDouble delta = mul_dd({q, 0.0}, r);

    return add_dd({q, 0.0}, delta);
}