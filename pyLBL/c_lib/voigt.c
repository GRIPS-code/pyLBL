#include <math.h>


void voigt(double * dwno, int start, int end, double nu, double alpha, double gamma,
           double sw, double * k)
{
    double rsqrpi = 1./sqrt(M_PI);
    double sqrln2 = sqrt(log(2.));
    double y0 = 1.5;
    double y0py0 = y0 + y0;
    double y0q = y0*y0;

    double repwid = sqrln2/alpha;
    double y = repwid*gamma;
    double yq = y*y;

    if (y >= 70.55)
    {
        /*region 0.*/
        int i;
        for (i=start; i<=end; ++i)
        {
            double xi = (dwno[i] - nu)*repwid;
            k[i] += sw*repwid*y/(M_PI*(xi*xi + yq));
        }
        return;
    }

    int rg1 = 1;
    int rg2 = 1;
    int rg3 = 1;

    double yrrtpi = y*rsqrpi; /*y/sqrt(pi)*/
    double xlim0 = sqrt(15100. + y*(40. - y*3.6)); /*y < 70.55*/
    double xlim1;
    if (y >= 8.425)
    {
        xlim1 = 0.;
    }
    else
    {
        xlim1 = sqrt(164. - y*(4.3 + y*1.8));
    }
    double xlim2 = 6.8 - y;
    double xlim3 = 2.4*y;
    double xlim4 = 18.1*y + 1.65;

    if (y <= 0.000001)
    {
        /*when y<10^-6*/
        xlim1 = xlim0; /*avoid w4 algorithm.*/
        xlim2 = xlim0;
    }

    double c[6] = {1.0117281, -0.75197147, 0.012557727,
                   0.010022008, -0.00024206814, 0.00000050084806};
    double s[6] = {1.393237, 0.23115241, -0.15535147,
                   0.0062183662, 0.000091908299, -0.00000062752596};
    double t[6] = {0.31424038, 0.94778839, 1.5976826,
                   2.2795071, 3.0206370, 3.8897249};
    double a0, d0, d2, e0, e2, e4, h0, h2, h4, h6;
    double p0, p2, p4, p6, p8, z0, z2, z4, z6, z8;
    double buf;
    double xp[6];
    double xm[6];
    double yp[6];
    double ym[6];
    double mq[6];
    double pq[6];
    double mf[6];
    double pf[6];

    int i;
    for (i=start; i<=end; ++i)
    {
        double xi = (dwno[i] - nu)*repwid; /*loop over all points*/
        double abx = fabs(xi); /*|x|*/
        double xq = abx*abx; /*x^2*/
        if (abx >= xlim0)
        {
            /*region 0 algorithm.*/
            buf = yrrtpi/(xq + yq);
        }
        else if (abx >= xlim1)
        {
            /*humlicek w4 region 1.*/
            if (rg1 != 0)
            {
                /*first point in region 1.*/
                rg1 = 0;
                a0 = yq + 0.5; /*region 1 y-dependents.*/
                d0 = a0 * a0;
                d2 = yq + yq - 1.;
            }
            double d = rsqrpi/(d0 + xq*(d2 + xq));
            buf = d*y*(a0 + xq);
        }
        else if (abx >= xlim2)
        {
            /*humlicek w4 region 2.*/
            if (rg2 != 0)
            {
                /*first point in region 2.*/
                rg2 = 0;
                h0 = 0.5625 + yq*(4.5 + yq*(10.5 + yq*(6.0 + yq))); /*region 2 y-deps.*/
                h2 = -4.5 + yq*(9.0 + yq*(6.0 + yq*4.0));
                h4 = 10.5 - yq*(6.0 - yq*6.0);
                h6 = -6.0 + yq* 4.0;
                e0 = 1.875 + yq*(8.25 + yq*(5.5 + yq));
                e2 = 5.25 + yq*(1.0 + yq*3.0);
                e4 = 0.75*h6;
            }
            double d = rsqrpi/(h0 + xq*(h2 + xq*(h4 + xq*(h6 + xq))));
            buf = d*y*(e0 + xq*(e2 + xq*(e4 + xq)));
        }
        else if (abx < xlim3)
        {
            /*humlicek w4 region 3*/
            if (rg3 != 0)
            {
                /*first point in region 3.*/
                rg3 = 0;
                z0 = 272.1014 + y*(1280.829 + y*(2802.870 + y*(3764.966 /*y-deps*/
                     + y*(3447.629 + y*(2256.981 + y*(1074.409 + y*(369.1989
                     + y*(88.26741 + y*(13.39880 + y)))))))));
                z2 = 211.678 + y*(902.3066 + y*(1758.336 + y*(2037.310
                     + y*(1549.675 + y*(793.4273 + y*(266.2987
                     + y*(53.59518 + y*5.0)))))));
                z4 = 78.86585 + y*(308.1852 + y*(497.3014 + y*(479.2576
                     + y*(269.2916 + y*(80.39278 + y*10.0)))));
                z6 = 22.03523 + y*(55.02933 + y*(92.75679 + y*(53.59518
                     + y*10.0)));
                z8 = 1.496460 + y*(13.39880 + y*5.0);
                p0 = 153.5168 + y*(549.3954 + y*(919.4955 + y*(946.8970
                     + y*(662.8097 + y*(328.2151 + y*(115.3772 + y*(27.93941
                     + y*(4.264678 + y*0.3183291))))))));
                p2 = -34.16955 + y*(-1.322256+ y*(124.5975 + y*(189.7730
                     + y*(139.4665 + y*(56.81652 + y*(12.79458
                     + y*1.2733163))))));
                p4 = 2.584042 + y*(10.46332 + y*(24.01655 + y*(29.81482
                     + y*(12.79568 + y*1.9099744))));
                p6 = -0.07272979  + y*(0.9377051+ y*(4.266322 + y*1.273316));
                p8 = 0.0005480304 + y*0.3183291;
            }
            double d = 1.7724538/(z0 + xq*(z2 + xq*(z4 + xq*(z6 + xq*(z8 + xq)))));
            buf = d*(p0 + xq*(p2 + xq*(p4 + xq*(p6 + xq*p8))));
        }
        else
        {
            /*humlicek cpf12 algorithm.*/
            double ypy0 = y + y0;
            double ypy0q = ypy0*ypy0;
            buf = 0.;
            int j;
            for (j=0; j<6; ++j)
            {
                double d = xi - t[j];
                mq[j] = d*d;
                mf[j] = 1./(mq[j] + ypy0q);
                xm[j] = mf[j]*d;
                ym[j] = mf[j]*ypy0;
                d = xi + t[j];
                pq[j] = d*d;
                pf[j] = 1./(pq[j] + ypy0q);
                xp[j] = pf[j]*d;
                yp[j] = pf[j]*ypy0;
            }
            if (abx <= xlim4)
            {
                /*humlicek cpf12 region i.*/
                for (j=0; j<6; ++j)
                {
                    buf += c[j]*(ym[j] + yp[j]) - s[j]*(xm[j] - xp[j]);
                }
            }
            else
            {
                /*humlicek cpf12 region ii.*/
                double yf = y + y0py0;
                for (j=0; j<6; ++j)
                {
                    buf += (c[j]*(mq[j]*mf[j] - y0*ym[j]) + s[j]*yf*xm[j])/(mq[j] + y0q)
                           + (c[j]*(pq[j]*pf[j] - y0*yp[j]) - s[j]*yf*xp[j])/(pq[j] + y0q);
                }
                buf = y*buf + exp(-xq);
            }
        }
        k[i] += sw*rsqrpi*repwid*buf;
    }
    return;
}
