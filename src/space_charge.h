#ifndef DEMOTRACK_SPACE_CHARGE_H__
#define DEMOTRACK_SPACE_CHARGE_H__

DEMOTRACK_STATIC DEMOTRACK_FN void cerrf(
    double const in_real, double const in_imag,
    double* out_real, double* out_imag );

DEMOTRACK_STATIC DEMOTRACK_FN void get_transv_field_gauss_round(
    double const sigma, double const Delta_x,
    double const Delta_y, double const x, double const y,
    double* Ex, double* Ey );

DEMOTRACK_STATIC DEMOTRACK_FN void get_transv_field_gauss_ellip(
    double const sigma_x, double const sigma_y,
    double const Delta_x, double const Delta_y,
    double const x, double const y,
    double* Ex, double* Ey );

DEMOTRACK_STATIC DEMOTRACK_FN void get_Ex_Ey_Gx_Gy_gauss(
    double const x, double const y,
    double const sigma_x, double const sigma_y,
    double const min_sigma_diff, int skip_Gs,
    double* Ex_ptr, double* Ey_ptr,
    double* Gx_ptr, double* Gy_ptr);

DEMOTRACK_STATIC DEMOTRACK_FN void SpaceChargeCoasting_track(
    BE_ARGPTR_DEC const SpaceChargeCoasting *const cavity,
    PARTICLE_ARGPTR_DEC Particle* p );

/* ************************************************************************* */

DEMOTRACK_INLINE void cerrf( double const in_real, double const in_imag,
    double* out_real, double* out_imag )
{
    /* This function calculates the double precision complex error fnct.
    based on the algorithm of the FORTRAN function written at CERN by K. Koelbig
    Program C335, 1970. See also M. Bassetti and G.A. Erskine, "Closed
    expression for the electric field of a two-dimensional Gaussian charge
    density", CERN-ISR-TH/80-06; */

    int n, nc, nu;
    double a_constant = ( double )1.12837916709551;
    double xLim = ( double )5.33;
    double yLim = ( double )4.29;
    double h, q, Saux, Sx, Sy, Tn, Tx, Ty, Wx, Wy, xh, xl, x, yh, y;
    double Rx [33];
    double Ry [33];

    x = fabs(in_real);
    y = fabs(in_imag);

    if (y < yLim && x < xLim){
        q = (1.0 - y / yLim) * sqrt(1.0 - (x / xLim) * (x / xLim));
        h  = 1.0 / (3.2 * q);
        nc = 7 + (int) (23.0 * q);
        xl = pow(h, (double) (1 - nc));
        xh = y + 0.5 / h;
        yh = x;
        nu = 10 + (int) (21.0 * q);
        Rx[nu] = 0.;
        Ry[nu] = 0.;
        for (n = nu; n > 0; n--){
            Tx = xh + n * Rx[n];
            Ty = yh - n * Ry[n];
            Tn = Tx*Tx + Ty*Ty;
            Rx[n-1] = 0.5 * Tx / Tn;
            Ry[n-1] = 0.5 * Ty / Tn;
            }
        /* .... */
        Sx = 0.;
        Sy = 0.;
        for (n = nc; n>0; n--){
            Saux = Sx + xl;
            Sx = Rx[n-1] * Saux - Ry[n-1] * Sy;
            Sy = Rx[n-1] * Sy + Ry[n-1] * Saux;
            xl = h * xl;
        };
        Wx = a_constant * Sx;
        Wy = a_constant * Sy;
    }
    else{
        xh = y;
        yh = x;
        Rx[0] = 0.;
        Ry[0] = 0.;
        for (n = 9; n>0; n--){
            Tx = xh + n * Rx[0];
            Ty = yh - n * Ry[0];
            Tn = Tx * Tx + Ty * Ty;
            Rx[0] = 0.5 * Tx / Tn;
            Ry[0] = 0.5 * Ty / Tn;
        };
        Wx = a_constant * Rx[0];
        Wy = a_constant * Ry[0];
    }
    if (y == 0.) {Wx = exp(-x * x);}
    if (in_imag < 0.){
        Wx =   2.0 * exp(y * y - x * x) * cos(2.0 * x * y) - Wx;
        Wy = - 2.0 * exp(y * y - x * x) * sin(2.0 * x * y) - Wy;
        if (in_real > 0.) {Wy = -Wy;}
    }
    else if (in_real < 0.) {Wy = -Wy;}

    *out_real = Wx;
    *out_imag = Wy;
}

DEMOTRACK_INLINE void get_transv_field_gauss_round(
    double const sigma, double const Delta_x,
    double const Delta_y, double const x, double const y,
    double* Ex_out, double* Ey_out )
{
  double r2, temp;
  double const PI = ( double )3.141592653589793;
  double const EPSILON_0 = ( double )8.854187817620e-12;

  r2 = (x-Delta_x)*(x-Delta_x)+(y-Delta_y)*(y-Delta_y);
  if (r2<1e-20) temp = sqrt(r2)/(2.*PI*EPSILON_0*sigma); //linearised
  else          temp = (1-exp(-0.5*r2/(sigma*sigma)))/(2.*PI*EPSILON_0*r2);

  (*Ex_out) = temp * (x-Delta_x);
  (*Ey_out) = temp * (y-Delta_y);
}

DEMOTRACK_INLINE void get_transv_field_gauss_ellip(
    double const sigma_x, double const sigma_y,
    double const Delta_x, double const Delta_y,
    double const x, double const y,
    double* Ex_out, double* Ey_out )
{
  double sigmax = sigma_x;
  double sigmay = sigma_y;

  double const abx = fabs(x - Delta_x);
  double const aby = fabs(y - Delta_y);

  double S, factBE, Ex, Ey;
  double etaBE_re, etaBE_im, zetaBE_re, zetaBE_im;
  double w_etaBE_re, w_etaBE_im, w_zetaBE_re, w_zetaBE_im;
  double expBE;

  double const SQRT_PI = ( double )1.7724538509055160273;
  double const EPSILON_0 = ( double )8.854187817620e-12;

  if (sigmax>sigmay){
    S = sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
    factBE = 1./(2.*EPSILON_0*SQRT_PI*S);

    etaBE_re = sigmay/sigmax*abx;
    etaBE_im = sigmax/sigmay*aby;

    zetaBE_re = abx;
    zetaBE_im = aby;

    cerrf(zetaBE_re/S, zetaBE_im/S , &(w_zetaBE_re), &(w_zetaBE_im));
    cerrf(etaBE_re/S, etaBE_im/S , &(w_etaBE_re), &(w_etaBE_im));

    expBE = exp(-abx*abx/(2*sigmax*sigmax)-aby*aby/(2*sigmay*sigmay));

    Ex = factBE*(w_zetaBE_im - w_etaBE_im*expBE);
    Ey = factBE*(w_zetaBE_re - w_etaBE_re*expBE);
  }
  else if (sigmax<sigmay){
    S = sqrt(2.*(sigmay*sigmay-sigmax*sigmax));
    factBE = 1./(2.*EPSILON_0*SQRT_PI*S);

    etaBE_re = sigmax/sigmay*aby;
    etaBE_im = sigmay/sigmax*abx;

    zetaBE_re = aby;
    zetaBE_im = abx;

    cerrf(zetaBE_re/S, zetaBE_im/S , &(w_zetaBE_re), &(w_zetaBE_im));
    cerrf(etaBE_re/S, etaBE_im/S , &(w_etaBE_re), &(w_etaBE_im));

    expBE = exp(-aby*aby/(2*sigmay*sigmay)-abx*abx/(2*sigmax*sigmax));

    Ey = factBE*(w_zetaBE_im - w_etaBE_im*expBE);
    Ex = factBE*(w_zetaBE_re - w_etaBE_re*expBE);
  }

  if((x - Delta_x)<0) Ex=-Ex;
  if((y - Delta_y)<0) Ey=-Ey;

  (*Ex_out) = Ex;
  (*Ey_out) = Ey;
}

DEMOTRACK_INLINE void get_Ex_Ey_Gx_Gy_gauss( double const x, double const y,
    double const sigma_x, double const sigma_y,
    double const min_sigma_diff, int skip_Gs,
    double* Ex_ptr, double* Ey_ptr,
    double* Gx_ptr, double* Gy_ptr)
{
    double const PI = ( double )3.141592653589793;
    double const EPSILON_0 = ( double )8.854187817620e-12;

    double Ex, Ey, Gx, Gy;

    if (fabs(sigma_x-sigma_y)< min_sigma_diff){

        double sigma = 0.5*(sigma_x+sigma_y);

        get_transv_field_gauss_round(sigma, 0., 0., x, y, &Ex, &Ey);

        if(skip_Gs){
          Gx = 0.;
          Gy = 0.;
        }
        else{
          Gx = 1/(2.*(x*x+y*y))*(y*Ey-x*Ex+1./(2*PI*EPSILON_0*sigma*sigma)
                            *x*x*exp(-(x*x+y*y)/(2.*sigma*sigma)));
          Gy = 1./(2*(x*x+y*y))*(x*Ex-y*Ey+1./(2*PI*EPSILON_0*sigma*sigma)
                            *y*y*exp(-(x*x+y*y)/(2.*sigma*sigma)));
        }

    }
    else{

        get_transv_field_gauss_ellip(
                sigma_x, sigma_y, 0., 0., x, y, &Ex, &Ey);

        double Sig_11 = sigma_x*sigma_x;
        double Sig_33 = sigma_y*sigma_y;

        if(skip_Gs){
          Gx = 0.;
          Gy = 0.;
        }
        else{
          Gx =-1./(2*(Sig_11-Sig_33))*(x*Ex+y*Ey+1./(2*PI*EPSILON_0)*\
                      (sigma_y/sigma_x*exp(-x*x/(2*Sig_11)-y*y/(2*Sig_33))-1.));
          Gy =1./(2*(Sig_11-Sig_33))*(x*Ex+y*Ey+1./(2*PI*EPSILON_0)*\
                      (sigma_x/sigma_y*exp(-x*x/(2*Sig_11)-y*y/(2*Sig_33))-1.));
        }
    }

    *Ex_ptr = Ex;
    *Ey_ptr = Ey;

    if( !skip_Gs )
    {
        *Gx_ptr = Gx;
        *Gy_ptr = Gy;
    }
}

DEMOTRACK_INLINE void SpaceChargeCoasting_track(
    BE_ARGPTR_DEC const SpaceChargeCoasting *const sc_elem,
    PARTICLE_ARGPTR_DEC Particle* p )
{
    double Ex;
    double Ey;
    double Gx;
    double Gy;

    double const charge = ( double )1.602176634e-19 * p->q0;
    double const p0c = ( double )1.602176634e-19 * p->p0c;

    double fact_kick = sc_elem->number_of_particles *
        sc_elem->length * p->chi * p->charge_ratio * charge * charge *
        ( ( double )1 - p->beta0 * p->beta0 );

    fact_kick /= sc_elem->circumference * p0c * p->beta0 * p->rvv;

    get_Ex_Ey_Gx_Gy_gauss( p->x - sc_elem->x_co, p->y - sc_elem->y_co,
        sc_elem->sigma_x, sc_elem->sigma_y, sc_elem->min_sigma_diff,
        ( int )1, &Ex, &Ey, &Gx, &Gy );

    p->px += fact_kick * Ex;
    p->py += fact_kick * Ey;
}

#endif /* DEMOTRACK_SPACE_CHARGE_H__ */
