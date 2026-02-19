#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dd_math.h"

#define STATE_DIM 4

typedef struct{
    DoubleDouble x, y, px, py;
} State;

DoubleDouble compute_energy(const State* s) {
    auto px2 = mul_dd(s->px, s->px);
    auto py2 = mul_dd(s->py, s->py);
    auto x2  = mul_dd(s->x,  s->x);
    auto y2  = mul_dd(s->y,  s->y);
    auto x2y2 = mul_dd(x2, y2);

    auto sum1 = add_dd(px2, py2);
    auto sum2 = add_dd(sum1, x2y2);

    return sum2;
}

void rhs(const State* s, DoubleDouble* dxdt) {
    dxdt[0] = s->px;                       
    dxdt[1] = s->py;                        
    dxdt[2] = mul_dd(mul_dd(s->x, s->y), -(s->y));
    dxdt[3] = mul_dd(mul_dd(s->y, s->x), -(s->x));
}

void rk4_step_dd(const State* y, DoubleDouble h, State* y_new) {

    DoubleDouble half_h   = mul_small(h, 0.5);
    DoubleDouble sixth_h  = mul_small(h, 1.0/6.0);
    DoubleDouble third_h  = mul_small(h, 1.0/3.0);

    

    DoubleDouble k1[4], k2[4], k3[4], k4[4];
    State temp;

    rhs(y, k1);

    temp.x  = add_dd(y->x,  mul_dd(k1[0], half_h));
    temp.y  = add_dd(y->y,  mul_dd(k1[1], half_h));
    temp.px = add_dd(y->px, mul_dd(k1[2], half_h));
    temp.py = add_dd(y->py, mul_dd(k1[3], half_h));
    rhs(&temp, k2);

    temp.x  = add_dd(y->x,  mul_dd(k2[0], half_h));
    temp.y  = add_dd(y->y,  mul_dd(k2[1], half_h));
    temp.px = add_dd(y->px, mul_dd(k2[2], half_h));
    temp.py = add_dd(y->py, mul_dd(k2[3], half_h));
    rhs(&temp, k3);

    temp.x  = add_dd(y->x, mul_dd(k3[0], half_h));
    temp.y  = add_dd(y->y,  mul_dd(k3[1], half_h));
    temp.px = add_dd(y->px, mul_dd(k3[2], half_h));
    temp.py = add_dd(y->py, mul_dd(k3[3], half_h));
    rhs(&temp, k4);


    for (int i = 0; i < STATE_DIM; ++i) {
        auto contrib1 = mul_dd(k1[i], sixth_h);
        auto contrib2 = mul_dd(k2[i], third_h);
        auto contrib3 = mul_dd(k3[i], third_h);
        auto contrib4 = mul_dd(k4[i], sixth_h);

        auto sum_contrib = add_dd(contrib1, add_dd(contrib2, add_dd(contrib3, contrib4)));

        if (i == 0)      y_new->x  = add_dd(y->x,  sum_contrib);
        else if (i == 1) y_new->y  = add_dd(y->y,  sum_contrib);
        else if (i == 2) y_new->px = add_dd(y->px, sum_contrib);
        else if (i == 3) y_new->py = add_dd(y->py, sum_contrib);
    }



}

int main()
{
    State y     = {};
    State y_new = {};

    double t_double = 0.0;
    double t_end    = 10.0;
    double h_double = 1e-8;

    DoubleDouble t     = to_dd(0.0);
    DoubleDouble h     = to_dd(h_double);
    DoubleDouble t_end_dd = to_dd(t_end);

    long long step = 0;

    y.x  = to_dd(1.0);
    y.y  = to_dd(0.0);
    y.px = to_dd(0.5);
    DoubleDouble px2     = mul_dd(y.px, y.px);
    DoubleDouble one     = to_dd(1.0);
    DoubleDouble arg_sqrt = sub_dd(one, px2);
    double arg_sqrt_approx = arg_sqrt.hi;
    y.py = to_dd(sqrt(arg_sqrt_approx));

    DoubleDouble E0 = compute_energy(&y);

    printf("Начальная энергия: %.15e   (lo = %+.3e)\n", E0.hi, E0.lo);

    while (true)
    {
        DoubleDouble diff = sub_dd(t_end_dd, t);
        if (diff.hi <= 0.0) break;

        step++;

        rk4_step_dd(&y, h, &y_new);

        DoubleDouble E_new   = compute_energy(&y_new);
        DoubleDouble delta_E = sub_dd(E_new, E0);
        double delta_E_abs   = fabs(delta_E.hi);

        if (delta_E_abs > 1e-11)
        {
            h = mul_dd(h, to_dd(0.5));
            continue;
        }
        else if (delta_E_abs < 1e-20)
        {
            h = mul_dd(h, to_dd(2.0));
        }

        E0 = E_new;

        if (step % 100'000'000 == 0)
        {
            printf("step %8lld   t ≈ %.6f   h = %.3e   E_err ≈ %.3e\n",
                   step, t.hi, h.hi, delta_E_abs);
        }

        DoubleDouble prod = mul_dd(y.y, y_new.y);
        if ((prod.hi <= 0.0 || (prod.hi == 0.0 && prod.lo < 0.0)) && t.hi > 1e-10)
        {
            DoubleDouble dy     = sub_dd(y_new.y, y.y);
            DoubleDouble frac   = mul_dd( sub_dd(to_dd(0.0), y.y), reciprocal_dd(dy)) ; 

            DoubleDouble neg_y   = -y.y; 
            DoubleDouble denom   = sub_dd(y_new.y, y.y);
            double frac_approx = neg_y.hi / denom.hi; 

            DoubleDouble t_event = add_dd(t, mul_dd(h, to_dd(frac_approx)));

            State y_event;
            y_event.x  = add_dd(y.x,  mul_dd( to_dd(frac_approx), sub_dd(y_new.x,  y.x)  ));
            y_event.y  = add_dd(y.y,  mul_dd( to_dd(frac_approx), sub_dd(y_new.y,  y.y)  ));
            y_event.px = add_dd(y.px, mul_dd( to_dd(frac_approx), sub_dd(y_new.px, y.px) ));
            y_event.py = add_dd(y.py, mul_dd( to_dd(frac_approx), sub_dd(y_new.py, y.py) ));

            printf("\n=== Пересечение y=0 обнаружено ===\n");
            printf("t_event ≈ %.10f\n", t_event.hi);
            printf("x       = %.10f\n", y_event.x.hi);
            printf("y       = %.10f  (lo = %+.2e)\n", y_event.y.hi, y_event.y.lo);
            printf("px      = %.10f\n", y_event.px.hi);
            printf("py      = %.10f\n", y_event.py.hi);
            printf("энергия = %.20e   (lo = %+.3e)\n", compute_energy(&y_event).hi,
                                                       compute_energy(&y_event).lo);

            break;
        }

        y = y_new;
        t = add_dd(t, h);
    }

    printf("\nИнтеграция завершена.\n");
    printf("Последнее t ≈ %.6f,  шагов = %lld,  последний h = %.3e\n",
           t.hi, step, h.hi);

    return 0;
}