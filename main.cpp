#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

const double ENERGY_ATOL = 1e-14;
const double ENERGY_RTOL = 1e-13;

#define STATE_DIM 4
#define MONO_DIM 16
#define FULL_DIM (STATE_DIM + MONO_DIM)

typedef struct
{
    double data[FULL_DIM];
} State;

double inline compute_energy(const State *s)
{
    return 0.5 * (s->data[2]*s->data[2] + s->data[3]*s->data[3] +
                  s->data[0]*s->data[0] * s->data[1]*s->data[1]);
}

void inline rhs(const State *s, double *dxdt)
{
    dxdt[0] = s->data[2];
    dxdt[1] = s->data[3];
    dxdt[2] = -s->data[0] * s->data[1] * s->data[1];
    dxdt[3] = -s->data[1] * s->data[0] * s->data[0];

    double xx = s->data[0];
    double yy = s->data[1];

    double J[4][4] = {
        {0,    0,    1,    0},
        {0,    0,    0,    1},
        {-yy*yy, -2*xx*yy, 0, 0},
        {-2*xx*yy, -xx*xx, 0, 0}
    };

    for (int col = 0; col < 4; col++)
    {
        int off = 4 + col * 4;
        double mx  = s->data[off + 0];
        double my  = s->data[off + 1];
        double mpx = s->data[off + 2];
        double mpy = s->data[off + 3];

        dxdt[off + 0] = J[0][0]*mx + J[0][1]*my + J[0][2]*mpx + J[0][3]*mpy;
        dxdt[off + 1] = J[1][0]*mx + J[1][1]*my + J[1][2]*mpx + J[1][3]*mpy;
        dxdt[off + 2] = J[2][0]*mx + J[2][1]*my + J[2][2]*mpx + J[2][3]*mpy;
        dxdt[off + 3] = J[3][0]*mx + J[3][1]*my + J[3][2]*mpx + J[3][3]*mpy;
    }
}

void inline rk4_step(const State *y, double h, State *y_new)
{
    double k1[FULL_DIM], k2[FULL_DIM], k3[FULL_DIM], k4[FULL_DIM];
    State temp;

    rhs(y, k1);
    for (int i = 0; i < FULL_DIM; i++)
        temp.data[i] = y->data[i] + 0.5 * h * k1[i];
    rhs(&temp, k2);

    for (int i = 0; i < FULL_DIM; i++)
        temp.data[i] = y->data[i] + 0.5 * h * k2[i];
    rhs(&temp, k3);

    for (int i = 0; i < FULL_DIM; i++)
        temp.data[i] = y->data[i] + h * k3[i];
    rhs(&temp, k4);

    for (int i = 0; i < FULL_DIM; i++)
        y_new->data[i] = y->data[i] + (h/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Использование: %s <файл_с_начальными_условиями>\n", argv[0]);
        return 1;
    }

    const char *input_filename = argv[1];
    FILE *f_in = fopen(input_filename, "r");
    if (!f_in) {
        fprintf(stderr, "Не удалось открыть %s\n", input_filename);
        return 1;
    }

    int n_cases = 0;
    if (fscanf(f_in, "%d", &n_cases) != 1 || n_cases <= 0) {
        fprintf(stderr, "Ошибка чтения количества наборов\n");
        fclose(f_in);
        return 1;
    }

    printf("Обнаружено %d наборов начальных условий\n", n_cases);

    for (int cas = 1; cas <= n_cases; cas++)
    {
        State y = {0}, y_new = {0};
        double t = 0.0;
        double t_end = 35.0 + 1e-7;
        double h = 1e-10;
        long long step = 0;
        bool near_zero = false;

        double x0, y0, px0, py0;
        if (fscanf(f_in, "%lf %lf %lf %lf", &x0, &y0, &px0, &py0) != 4) {
            fprintf(stderr, "Ошибка чтения набора %d\n", cas);
            break;
        }

        y.data[0] = x0; y.data[1] = y0;
        y.data[2] = px0; y.data[3] = py0;

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                y.data[4 + j*4 + i] = (i == j) ? 1.0 : 0.0;

        double E0 = compute_energy(&y);
        printf("\n[%2d] Начало:  x=%.12e  y=%.12e  px=%.12e  py=%.12e   E=%.16e\n",
               cas, x0, y0, px0, py0, E0);

        char filename[64];
        snprintf(filename, sizeof(filename), "trajectory_%d.txt", cas);

        FILE *fp = fopen(filename, "w");
        if (!fp) {
            fprintf(stderr, "Не удалось создать %s\n", filename);
            continue;
        }

        fprintf(fp, "# t                  x                    y                    px                   py                   energy\n");
        fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.16e\n",
                t, y.data[0], y.data[1], y.data[2], y.data[3], E0);

        double step_factor;
        uint64_t istep = 0x3FF0010000000000ULL;
        memcpy(&step_factor, &istep, sizeof(double));
        double step_factor_inv = 1.0 / step_factor;

        while (t < t_end)
        {
            rk4_step(&y, h, &y_new);

            double E_new = compute_energy(&y_new);
            double delta_E = fabs(E_new - E0);

            if (delta_E > 1e-13) {
                h *= step_factor_inv;
                continue;
            }
            else if (delta_E < 1e-17 && !near_zero) {
                h *= step_factor;
                continue;
            }

            E0 = E_new;

            if (step % 100 == 0 && step > 0) {
                fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.16e\n",
                    t, y.data[0], y.data[1], y.data[2], y.data[3], E_new);
            }

            bool crossed_y  = (y.data[1] * y_new.data[1] <= 0.0) && t > 1e-10;
            bool crossed_px = (y.data[2] * y_new.data[2] <= 0.0) && t > 1e-10;

            if (crossed_y || crossed_px)
            {
                if ((crossed_y  && fabs(y_new.data[1]) > 1e-100) ||
                    (crossed_px && fabs(y_new.data[2]) > 1e-100))
                {
                    near_zero = true;
                    h *= step_factor_inv;
                    continue;
                }

                const char *event = crossed_y ? "y=0" : "px=0";

                printf("\n[%2d] === Пересечение %s ===\n", cas, event);
                printf("   t ≈ %.19f\n", t + h);
                printf("   x  = %.19f\n", y_new.data[0]);
                printf("   y  = %.19f\n", y_new.data[1]);
                printf("   px = %.19f\n", y_new.data[2]);
                printf("   py = %.19f\n", y_new.data[3]);
                printf("   E  = %.19e   ΔE = %.2e\n", E_new, delta_E);

                printf("   Монодромия:\n");
                for (int r = 0; r < 4; r++) {
                    printf("   [ ");
                    for (int c = 0; c < 4; c++) {
                        printf("%12.6e", y.data[4 + c*4 + r]);
                        if (c < 3) printf(", ");
                    }
                    printf("]");
                    if (r < 3) printf(",\n");
                    else printf("\n");
                }

                fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.16e\n",
                        t + h, y_new.data[0], y_new.data[1], y_new.data[2], y_new.data[3], E_new);

                near_zero = false;
            }

            y = y_new;
            t += h;
            step++;
        }

        fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.16e   # конец интеграции\n",
                t, y.data[0], y.data[1], y.data[2], y.data[3], E0);

        fclose(fp);

        printf("[%2d] Завершено. t=%.6f, шагов=%lld, файл: %s\n",
               cas, t, step, filename);
    }

    fclose(f_in);
    printf("Все наборы обработаны.\n");

    return 0;
}