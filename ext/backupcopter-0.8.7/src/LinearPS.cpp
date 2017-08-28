#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "GrowthFunction.h"
#include "LinearPS.h"


LinearPS::LinearPS(const Cosmology& C_, real z_)
    : C(C_), z(z_), Pnw(C_, z_, EisensteinHu)
{
    GrowthFunction D(C);
    real A = pow2(C.delta_H) * pow2(D(z)/D(0));

    real T0 = C.Ti[0];
    int N = C.ki.size();

    /* Calculate P(k) */
    array p(N);
    for(int i = 0; i < N; i++) {
        real k = C.ki[i];
        real T = C.Ti[i];
        p[i] = 2*M_PI*M_PI/pow3(k) * A * pow(2997.925*k, 3+C.n) * pow2(T/T0);
    }
    pk = LinearSpline(C.ki, p);

    k0 = C.ki[0];
    k1 = C.ki[N-1];
    p0 = p[0];
    p1 = p[N-1];

    if(k0 <= 0.)
        error("LinearPS: k0 = %g\n", k0);
}

real LinearPS::Evaluate(real k) const {
    if(k <= 0)
        return 0;
    else if(k <= k0)
        return p0 * Pnw(k)/Pnw(k0);     // Eisenstein-Hu, scaled to match splined P(k) at k = k0
    else if(k >= k1)
        return p1 * Pnw(k)/Pnw(k1);     // Eisenstein-Hu, scaled to match splined P(k) at k = k1
    else
        return pk(k);
}
