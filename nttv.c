/*  My attempt at writing N-body integrator code 
    Based on legacy code written by Dan Fabrycky

    compiles with: gcc -O3 -o mynttv -lgsl -lgslcblas nttvfb.c

*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>



/* INTEGRATION PARAMETERS */
#define NPL 4                       /* Number of planets to integrate */
#define DY 1e-13                    /* Error allowed in parameter values per timestep */
#define HStart 1e-5                 /* Timestep to start; if you get NaNs, try reducing this */
                                    /* Change these to be arguments to the main function */

/* PHYSICAL CONSTANTS */
#define G 2.9591220363e-4           /* Newton's constant, in AU^3/days^2 */
                                    /* 2.9591439153e-4 (Laughlin) | 2.95912200e-4 (Dan-Eric) */


/* define operation to set vectors */
int seteq (double ynew[], const double y[])
{
  int i;
  for(i=0;i<6*NPL;i++) ynew[i]=y[i];
  return 0;
}


/* define the function */
int func (double t, const double y[], double f[], void *params)
{
    printf("running FUNCTION");

    int i,j;
    int i1,i2;
    double *masses = (double *)params;
    double gmc1, gmc2, gm1, gm2;
    double rc1m3, rc2m3, r12m3;

    /* Interaction between each planet and central star */
    for(i=0; i<NPL; i++) {
        gmc1  = G*(masses[0]+masses[i+1]);                                      /* G*(Mstar+Mpl)    */
        rc1m3 = pow(pow(y[i*6+0],2)+pow(y[i*6+1],2)+pow(y[i*6+2],2),-3.0/2);    /* rpl^-3           */

        for(j=0; j<3; j++) {
            f[i*6+j] = y[i*6+3+j];                                              /* x dot = v        */
            f[i*6+3+j] = -gmc1*y[i*6+j]*rc1m3;                                  /* Pull of the star */
        }   
    }
 
    /* Interaction between each pair of planets. */
    for(i1=0; i1<NPL-1; i1++) {
        gm1   = G*masses[i1+1];
        rc1m3 = pow(pow(y[i1*6+0],2)+pow(y[i1*6+1],2)+pow(y[i1*6+2],2),-3.0/2);

        for(i2=i1+1; i2<NPL; i2++) {
            gm2   = G*masses[i2+1];
            rc2m3 = pow(pow(y[i2*6+0],2)+pow(y[i2*6+1],2)+pow(y[i2*6+2],2),-3.0/2);
            r12m3 = pow(pow(y[i1*6+0]-y[i2*6+0],2)+pow(y[i1*6+1]-y[i2*6+1],2)+pow(y[i1*6+2]-y[i2*6+2],2),-3.0/2);
  
            for(j=0; j<3; j++) f[i1*6+3+j] += -gm2*( (y[i1*6+j]-y[i2*6+j])*r12m3 + y[i2*6+j]*rc2m3 );
            for(j=0; j<3; j++) f[i2*6+3+j] += -gm1*( (y[i2*6+j]-y[i1*6+j])*r12m3 + y[i1*6+j]*rc1m3 );
        }
    }

    return GSL_SUCCESS;
}



/* define the jacobian */
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    printf("running JACOBIAN");
    int i,j;
    double * masses = (double *)params;
    double gmc1 = G*(masses[0]+masses[1]);
    double gmc2 = G*(masses[0]+masses[2]);
    double gm1  = G*masses[1];
    double gm2  = G*masses[2];

    double rc1m3 = pow( pow(y[0],2)      + pow(y[1],2)      + pow(y[2],2),      -3.0/2);
    double rc2m3 = pow( pow(y[6],2)      + pow(y[7],2)      + pow(y[8],2),      -3.0/2);
    double r12m3 = pow( pow(y[0]-y[6],2) + pow(y[1]-y[7],2) + pow(y[2]-y[8],2), -3.0/2);
    double rc1m5 = pow( pow(y[0],2)      + pow(y[1],2)      + pow(y[2],2),      -5.0/2);
    double rc2m5 = pow( pow(y[6],2)      + pow(y[7],2)      + pow(y[8],2),      -5.0/2);
    double r12m5 = pow( pow(y[0]-y[6],2) + pow(y[1]-y[7],2) + pow(y[2]-y[8],2), -5.0/2);

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 6*NPL, 6*NPL);
    gsl_matrix *m            = &dfdy_mat.matrix; 

    /* Planet 1. */
    for(i=0;i<=2;i++){
        for(j=0;j<=2;j++)  gsl_matrix_set (m,i,j, 0.0);
        for(j=3;j<=5;j++)  gsl_matrix_set (m,i,j, 1.0);
        for(j=6;j<=11;j++) gsl_matrix_set (m,i,j, 0.0);
    }
    for(i=3;i<=5;i++){
        for(j=0;j<=2;j++){
            if(i==j+3) {
	            gsl_matrix_set 
	            (m,i,j, -gmc1*rc1m3 + 3*gmc1*y[j]*y[j]*rc1m5 - gmc2*(r12m3-3*y[j]*(y[j]-y[j+6])*r12m5));
            } else {
	            gsl_matrix_set 
	            (m,i,j, 3*gmc1*y[j]*y[i-3]*rc1m5 + gmc2*3*y[j]*(y[i-3]-y[i+3])*r12m5);
            }
        }
        for(j=3;j<=5;j++) gsl_matrix_set (m,i,j, 0.0);
        for(j=6;j<=8;j++){
            if(i==j-3) {
	            gsl_matrix_set 
	            (m,i,j, -gmc2*(-r12m3 + 3*y[j]*(y[j-6]-y[j])*r12m5 + rc2m3 - 3*y[j]*y[j]*rc2m5));
            } else {
	            gsl_matrix_set 
	            (m,i,j, -gmc2*(3*y[j]*(y[i-3]-y[i+3])*r12m5) - 3*y[j]*y[i+3]*rc2m5);
            }
        }
        for(j=9;j<=11;j++) gsl_matrix_set (m,i,j, 0.0);
    }

    /* Planet 2. */
    for(i=6;i<=8;i++){
        for(j=0;j<=8;j++)  gsl_matrix_set (m,i,j, 0.0);
        for(j=9;j<=11;j++) gsl_matrix_set (m,i,j, 1.0);
    }
    for(i=9;i<=11;i++){
        for(j=0;j<=2;j++){
            if(i==j-9) {
	            gsl_matrix_set 
	            (m,i,j, -gmc1*(-r12m3 + 3*y[j]*(y[j+6]-y[j])*r12m5 + rc1m3 - 3*y[j]*y[j]*rc1m5));
            } 
            else {
	            gsl_matrix_set 
	            (m,i,j, -gmc1*(3*y[j]*(y[i-3]-y[i-9])*r12m5) - 3*y[j]*y[i-9]*rc1m5);
            }
        }
        for(j=3;j<=5;j++) gsl_matrix_set (m, i, j, 0.0);
        for(j=6;j<=8;j++){
            if(i==j+3) {
	            gsl_matrix_set 
	            (m,i,j, -gmc2*rc2m3 + 3*gmc2*y[j]*y[j]*rc2m5 - gmc1*(r12m3 - 3*y[j]*(y[j]-y[j-6])*r12m5));
            }
            else {
	            gsl_matrix_set 
	            (m,i,j, 3*gmc2*y[j]*y[i-3]*rc2m5 + gmc1*3*y[j]*(y[i-3]-y[i-9])*r12m5);
            }
        }
        for(j=9;j<=11;j++) gsl_matrix_set (m,i,j, 0.0);
    }
    for(i=0;i<=11;i++) dfdt[i] = 0.0;

    return GSL_SUCCESS;
}


/* now write the integrator */
double *integrator(double *state, double *mu, double *tbounds)
{
    printf("running INTEGRATOR");

    /* inputs   *xyz     = (NPL x 6) state vector for planets [x,y,z,xdot,ydot,zdot...]
                *mu      = (NPL + 1) vector of dynamical masses, starting with central star
                *tbounds = (3) vector of times [t0,t1,tepoch]
    */
    
    /* set up integrator information */
    const gsl_odeiv_step_type   *T = gsl_odeiv_step_rk8pd;
          gsl_odeiv_step        *s = gsl_odeiv_step_alloc(T,6*NPL);
          gsl_odeiv_control     *c = gsl_odeiv_control_y_new(DY,0.0);
          gsl_odeiv_evolve      *e = gsl_odeiv_evolve_alloc(6*NPL);
    
    /* set up system */      
    gsl_odeiv_system sys = {func, jac, 6*NPL, mu};

    /* define integration bounds */
    double t;                       
    double t0=tbounds[0];           /* backward time to end integration; typically just before first RV or TTV datum */
    double t1=tbounds[1];           /* forward time to end integration; typically just after first RV or TTV datum */
    double tepoch=tbounds[2];       /* DCF: "tepoch" --> "tstart" */

    /* define integration quantities */
    long i;
    
    double y[6*NPL], yin[6*NPL];
    for(i=0;i<(6*NPL);i++) yin[i] = state[i];

    double mtot=mu[0];        
    for(i=0;i<NPL;i++) mtot += mu[i+1];
    
    double yhone[6*NPL], dydt[6*NPL], yerrhone[6*NPL];
    double h,thone,tstep;

    int transitcount[NPL], pl;
    double dist,vtrans;
    double tbv[4*NPL];                          /* transit time, impact parameter, velocity */
    double dp[NPL], ddpdt, dpold[NPL];          /* dp = dot-product = x(dot)v */
    double ps2, dist2;                          /* ps2 = projected separation squared; dist2 = distance squared */

    /* set up forward integration */
    seteq(y,yin);
    t = tepoch;
    h = HStart;
    for(pl=0;pl<NPL;pl++) {
        transitcount[pl] = -1;
        dp[pl]=y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];
    }
    
    /* do forward integration */
    while(t < t1){
 
        int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, y);

        if (status != GSL_SUCCESS)
            break;

        dist2 = pow(y[0],2)+pow(y[1],2)+pow(y[2],2);
    
        /* cyle through the planets, searching for a transit */
        for(pl=0; pl<NPL; pl++) {
            dpold[pl] = dp[pl];
	        dp[pl]    = y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];
	        ps2       = y[0+pl*6]*y[0+pl*6]+y[1+pl*6]*y[1+pl*6];

	        /* ps2 < dist2/5 --> on face of star | y[0] --> primary eclipse */	        
            if(dp[pl]*dpold[pl] <= 0 && ps2 < dist2/5 && y[2+pl*6] < 0) { 
	            transitcount[pl] += 1;
	            seteq(yhone, y);
	            thone = t;
	            
	            /* WHY IS THERE A FOR LOOP HERE? */
	            for(i=0; i<5; i++) {
	                func (thone, yhone, dydt, mu);
	                ddpdt = dydt[0+pl*6]*dydt[0+pl*6] + dydt[1+pl*6]*dydt[1+pl*6] 
	                         + y[0+pl*6]*dydt[3+pl*6] +    y[1+pl*6]*dydt[4+pl*6];
	                tstep = - dp[pl] / ddpdt;
	                gsl_odeiv_step_apply (s, thone, tstep, yhone, yerrhone, dydt, NULL, &sys);
	                thone += tstep;
	                dp[pl]=yhone[0+pl*6]*yhone[3+pl*6]+yhone[1+pl*6]*yhone[4+pl*6];
	            }
	            
	            dp[pl] = y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];
	            dist   = sqrt(pow(yhone[0+pl*6],2)+pow(yhone[1+pl*6],2));
	            vtrans = sqrt(pow(yhone[3+pl*6],2)+pow(yhone[4+pl*6],2));

                tbv[0+pl*4] = transitcount[pl];     /* transit number */
                tbv[1+pl*4] = thone;                /* transit time */
                tbv[2+pl*4] = dist;                 /* distance from star; ~impact parameter */
                tbv[3+pl*4] = vtrans;               /* velocity at transit */
            }
        }
    }

    /* setup backward integration */
    seteq(y, yin);
    t = tepoch;
    h = -HStart;
    for(pl=0;pl<NPL;pl++) {
        transitcount[pl] = 0;
        dp[pl]=y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];
    }

    /* do forward integration */
    while (t > t0+1) {      

        int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t0, &h, y);

        if (status != GSL_SUCCESS)
            break;

        dist2 = pow(y[0],2)+pow(y[1],2)+pow(y[2],2);
    
        /* Cycle through the planets, searching for a transit. */
        for(pl=0; pl<NPL; pl++) {  
	        dpold[pl] = dp[pl];
	        dp[pl]    = y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];
	        ps2       = y[0+pl*6]*y[0+pl*6]+y[1+pl*6]*y[1+pl*6];
	        
	        /* ps2 < dist2/5 --> on face of star | y[0] --> primary eclipse */
            if(dp[pl]*dpold[pl] <= 0 && ps2 < dist2/5 && y[2+pl*6] < 0) { 
	            transitcount[pl] += 1;
	            seteq(yhone, y);
	            thone = t;	        

	            /* WHY IS THERE A FOR LOOP HERE? */
	            for(i=0; i<5; i++) {
                    func (thone, yhone, dydt, mu);
                    ddpdt = dydt[0+pl*6]*dydt[0+pl*6] + dydt[1+pl*6]*dydt[1+pl*6]
                           + y[0+pl*6]*dydt[3+pl*6] +    y[1+pl*6]*dydt[4+pl*6];
                    tstep = - dp[pl] / ddpdt;
                    gsl_odeiv_step_apply (s, thone, tstep, yhone, yerrhone, dydt, NULL, &sys);
                    thone += tstep;
                    dp[pl]=yhone[0+pl*6]*yhone[3+pl*6]+yhone[1+pl*6]*yhone[4+pl*6];
	            }

                dp[pl] = y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];
                dist   = sqrt(pow(yhone[0+pl*6],2)+pow(yhone[1+pl*6],2));
                vtrans = sqrt(pow(yhone[3+pl*6],2)+pow(yhone[4+pl*6],2));
                
                tbv[0+pl*4] = transitcount[pl];     /* transit number */
                tbv[1+pl*4] = thone;                /* transit time */
                tbv[2+pl*4] = dist;                 /* distance from star; ~impact parameter */
                tbv[3+pl*4] = vtrans;               /* velocity at transit */
            }
        }
    }                

    /* free ode parameters */
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
	
    /* dummy array to return -- I'll deal with returning the full array later */
    double *ret = malloc(6 * sizeof(double));
    
    ret[0] = 1;
    ret[1] = 2;
    ret[2] = 3;
    ret[3] = 5;
    ret[4] = 8;
    ret[5] = 21;

    return ret;
}


/* function to free a pointer */
void myfree(void *pointer)
{
    free(pointer);
}
    
