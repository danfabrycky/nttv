/*  Code to solve Van der Pol oscillator equation

    2nd order:  u''(t) + mu * u'(t) * (u(t)^2 - 1) + u(t) = 0   

    1st order:  u' = v
                v' = -u + mu * v * (1 - u^2)

    Define the n-diminsional 1st order ODE as dy/dt = f(t,y)
        y(t) = [y1(t), y2(t) ... yn(t)]
        f(t) = [f1(t), f2(t) ... fn(t)]
        
    In particular here y(t) = [u(t),v(t)] and f(t) = [u'(t),v'(t)]
*/


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>


/* define the function */
int func (double t, const double y[], double f[], void *params)
{
    double mu = *(double *)params;
    f[0] = y[1];
    f[1] = -y[0] + mu*y[1]*(1.0-y[0]*y[0]);

    return GSL_SUCCESS;
}


/* define the jacobian */
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    double mu = *(double *)params;
    
    /* declare matrix */
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy,2,2);
    gsl_matrix *m = &dfdy_mat.matrix;
    
    /* set up matrix */
    gsl_matrix_set (m,0,0, 0.0);
    gsl_matrix_set (m,0,1, 1.0);
    gsl_matrix_set (m,1,0, -1.0-2.0*mu*y[0]*y[1]);
    gsl_matrix_set (m,1,1, mu*(1.0-y[0]*y[0]));
    
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    
    return GSL_SUCCESS;
}


/* now write the solver */
double *solver(void)
{
    /* set up step information */
    const gsl_odeiv_step_type   *T = gsl_odeiv_step_rk8pd;
          gsl_odeiv_step        *s = gsl_odeiv_step_alloc(T,2);
          gsl_odeiv_control     *c = gsl_odeiv_control_y_new(1e-6,0.0);
          gsl_odeiv_evolve      *e = gsl_odeiv_evolve_alloc(2);
    
    /* set up system */
    double mu = 10;      
    gsl_odeiv_system sys = {func, jac, 2, &mu};
    
    double t0 = 0.0, t1 = 100.0;
    double t  = t0;
    double y0[2] = {1.0,0.0};
    double y[2]  = {y0[0],y0[1]};
    double h = 1e-6;
    
    /* printf("start:\t%.5e, %.5e, %.5e\n", t, y[0], y[1]); */
    
    /* now run the ODE solver */
    while (t < t1) {
        int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, y);
        
        if (status != GSL_SUCCESS)
            break;
    }

    /* printf("end:\t%.5e, %.5e, %.5e\n", t, y[0], y[1]); */
   
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);   

    double *ret = malloc(6 * sizeof(double));
    
    ret[0] = t0;
    ret[1] = y0[0];
    ret[2] = y0[1];
    ret[3] = t;
    ret[4] = y[0];
    ret[5] = y[1];

/*    double ret[6] = {t0,y0[0],y0[1],t,y[0],y[1]};  */
    return ret;

}


/* function to free a pointer */
void myfree(void *pointer)
{
    free(pointer);
}



/* here's the main program */
/*
int main(void)
{
    double *r = solver();
    
    printf("%u\n", r);
    free(r);
        
    return 0;
}
*/
    
