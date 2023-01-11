
#include <math.h>



#include "def_usrmod.h"
#include "def_ind.h"

#define  NMOS   1
#define  NP     0
#define  NRC    0
#define  NRCE   0

#define  NXD    4
#define  NXA    0
#define  NU     2
#define  NPR    0

#define  k10    1.287e12
#define  k20    1.287e12
#define  k30    9.043e09
#define  E1     (-9758.3)
#define  E2     (-9758.3)
#define  E3     (-8560)
#define  H1     4.2
#define  H2     (-11.0)
#define  H3     (-41.85)
#define  rho    0.9342
#define  Cp     3.01
#define  kw     4032
#define  AR     0.215
#define  VR     10
#define  mK     5.0
#define  CPK    2.0

#define  cA0      5.1
#define  theta0   104.9

#define  FFs      14.19  /* Feed Flow normed by VR: (dV/dt  / VR)*/
#define  QdotKs   (-1113.5)

#define  cAs      2.1402105301746182E+00
#define  cBs      1.0903043613077321E+00
#define  thetas   1.1419108442079495E+02
#define  thetaKs  1.1290659291045561E+02


/*
2.14  
1.09  
114.2 
112.9 
*/

#define  TIMEUNITS_PER_HOUR 3600.0

#define  Q11    0.2
#define  Q22    1.0
#define  Q33    0.5
#define  Q44    0.2

#define  R11      0.5000
#define  R22      0.0000005


#define  P11    3278.78    
#define  P21    1677.31
#define  P31    681.02
#define  P41    271.50

#define  P12    1677.31
#define  P22    919.78
#define  P32    344.19
#define  P42    137.27

#define  P13    681.02
#define  P23    344.19
#define  P33    172.45
#define  P43    65.53

#define  P14    271.50
#define  P24    137.27
#define  P34    65.53
#define  P44    29.28


#define  R_OMEGA  90.0



static void ffcn(double *t, double *xd, double *xa, double *u, 
  double *p, double *rhs, double *rwh, long *iwh, long *info)
{ 
double k1, k2, k3;

k1=k10*exp(E1/(273.15 +xd[2]));
k2=k20*exp(E2/(273.15 +xd[2]));
k3=k30*exp(E3/(273.15 +xd[2]));

rhs[0] =(1/TIMEUNITS_PER_HOUR)*
(u[0]*(cA0-xd[0]) - k1*xd[0] - k3*xd[0]*xd[0]); 

rhs[1] =(1/TIMEUNITS_PER_HOUR)* 
(- u[0]*xd[1] + k1*xd[0] - k2*xd[1]); 

rhs[2] =(1/TIMEUNITS_PER_HOUR)*
(u[0]*(theta0-xd[2]) - (1/(rho*Cp)) *(k1*xd[0]*H1 + k2*xd[1]*H2 + k3*xd[0]*xd[0]*H3)
+(kw*AR/(rho*Cp*VR))*(xd[3] -xd[2])); 

rhs[3] =(1/TIMEUNITS_PER_HOUR)*
((1/(mK*CPK))*(u[1] + kw*AR*(xd[2]-xd[3])));


}



static void rdfcn_s(double *ts, double *sd, double *sa, double *u, 
  double *p, double *pr, double *res, long *dpnd, long *info)
{
  double alpha=0.5;
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, 0, 0, 0);
    return;
  }
  
  res[0] = sd[0]-(0.0*alpha+(1-alpha)*cAs);
  res[1] = sd[1]-(0.0*alpha+(1-alpha)*cBs);
  res[2] = sd[2]-(85.0*alpha+(1-alpha)*thetas);
  res[3] = sd[3]-(85.0*alpha+(1-alpha)*thetaKs);

}

static void rdfcn_e(double *ts, double *sd, double *sa, double *u, 
  double *p, double *pr, double *res, long *dpnd, long *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, 0, 0, 0);
    return;
  }

  res[0] = (sd[0]-cAs);
  res[1] = (sd[1]-cBs);
  res[2] = (sd[2]-thetas);
  res[3] = (sd[3]-thetaKs);


}
static void lsqfcn(double *ts, double *sd, double *sa, double *u, 
  double *p, double *pr, double *res, long *dpnd, long *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0);
    return;
  }

  res[0] = (sd[0]-cAs)*sqrt(Q11);
  res[1] = (sd[1]-cBs)*sqrt(Q22);
  res[2] = (sd[2]-thetas)*sqrt(Q33);
  res[3] = (sd[3]-thetaKs)*sqrt(Q44);

  res[4] = (u[0]-FFs)*sqrt(R11);
  res[5] = (u[1]-QdotKs)*sqrt(R22);

}

static void lfcn(double *t, double *xd, double *xa, double *u, 
  double *p, double *lval, double * rwh, long *iwh, long *info)
{
  static double res[6];
  static double pr[1];
  static long dpnd[1];
  long i;

  lsqfcn(t, xd, xa, u, 
	 p, pr, res, dpnd, info);

  *lval=0;
  for(i=0;i<6;i++)
    *lval+=res[i]*res[i];
}

void def_model(void)
{
  def_mdims(NMOS, NP, NRC, NRCE);


  /*/

  def_msolver(1, def_RKF45S);

  def_mstage(
    0, 
    NXD, NXA, NU, 
    NULL, lfcn, 
    0, 0, 0, NULL, ffcn, NULL, 
    NULL, NULL, 
    0
    );

  /*/
  def_msolver(1, def_adfDAESOL);

  def_mstage(
    0, 
    NXD, NXA, NU, 
    NULL, NULL, 
    0, 0, 0, NULL, ffcn, NULL, 
    NULL, NULL, 
    0
    );

  def_lsq(0, "c", NXD+NU, lsqfcn);





/*    def_lsq(1, "c", NXD+NU, lsqfcn); */
/*       def_mpc(0, "Start Point", NPR, NXD, NXD, rdfcn_s, NULL); */

}



