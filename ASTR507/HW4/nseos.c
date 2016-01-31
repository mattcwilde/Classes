/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **/
/*    ANALYTICAL REPRESENTATIONS OF NEUTRON-STAR EQUATIONS OF STATE     */
/*      Remarks and suggestions are welcome. Please send them to        */
/*         Alexander Potekhin <palex@astro.ioffe.ru>                    */
/*   For theoretical background and references see                      */
/*   P.Haensel & A.Y. Potekhin, Astron. Astrophys., 428, 191 (2004)     */
/*       and http://www.ioffe.ru/astro/NSG/NSEOS/                       */
/*   Last update: 31.08.04                                              */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **/

/* NB: You may need to compile with "-lm" to include "math.h" library.
       For example: cc nseos.c -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Prototypes of functions: */
int NSEOSFIT(int MODE, int KEOS, double *XN, double *RHO, double *H1,
 double *P, double *Gamma);
double FitRofN(int KEOS, double XN),
       FitNofR(int KEOS, double RHO),
       FitHEOS(int KEOS, double H1),
       FERMI(double X),
       FitGammaP(int KEOS, double RLG, double *PLG);

/*   ----------------------   MAIN block   ---------------------------  */
/*       This is auxiliary MAIN program for input/output purposes.      */
/*                You can change it or write your own.                  */
/*  Alternatively, you can just delete this MAIN block and link the set */
/*              of the remaining functions with yours.                  */
/*     Calculations are performed in the function NSEOSFIT below.       */
/*   -----------------------------------------------------------------  */

int main(void)
{
int MODE, KEOS, IRUN;
double XN, RHO, H1, P, Gamma;
   printf(" Arguments:\n");
   printf(" (1) number density of baryons n (in fm^{-3}),\n");
   printf(" (2) mass density of baryons rho (in g/cc),\n");
   printf(" (3) excess enthalpy per baryon h1=h/h0-1, where h0=m_0*c^2,\n");
   printf(" (4) pressure P (in dyn/cm^2),\n");
   printf(" (5) adiabatic exponent Gamma.\n");
   printf(" For one input argument (1, 2, or 3), other arguments are fitted.\n");
L0:  printf(" Input mode: 1 (n), 2 (rho), or 3 (h1) ? ");
   scanf("%d",&MODE);
   if (MODE < 1 || MODE > 3) {goto L0;};
L1:  printf(" EOS: 1 (SLy) or 2 (FPS) [enter 0 to stop] ? ");
   scanf("%d",&KEOS);
   if (KEOS == 0) {exit(0);};
   if (KEOS < 0 || KEOS > 2) {goto L1;};
   printf(" mode EOS#   n         rho        h1         P         Gamma\n");
     IRUN=0;
L10: IRUN=IRUN+1;
   if (MODE==1) printf(" n");
   else if (MODE==2) printf(" rho");
   else if (MODE==3) printf(" h1");
   if (IRUN==1) printf("  ( < 0 to terminate)");
   printf(": ");
   if (MODE==1) {
      scanf("%lf",&XN); 
     if (XN<=0.) goto L1;
   } else if (MODE==2) {
      scanf("%lf",&RHO); 
     if (RHO<=0.) goto L1;
   } else if (MODE==3) {
      scanf("%lf",&H1); 
     if (H1<=0.) goto L1;
   } else {
      printf(" Unknown MODE\n"); exit(0);
   };
   NSEOSFIT(MODE,KEOS,&XN,&RHO,&H1,&P,&Gamma);
   printf("%4d%4d%11.3E%11.3E%11.3E%11.3E%11.3E\n",MODE,KEOS,XN,RHO,H1,P,Gamma);
   goto L10;
}
/******************** END OF main ******************/


int NSEOSFIT(int MODE, int KEOS, double *XN, double *RHO, double *H1,
 double *P, double *Gamma)
/*                                                       Version 31.08.04
 Neutron-star EOS fitting.
 Arguments:
 MODE regulates input/output (see below),
 KEOS=1 for SLy EOS, 2 for FPS EOS,
 XN  - number density of baryons n (in fm^{-3}),
 RHO - mass density of baryons rho (in g/cc),
 H1  - excess enthalpy per baryon h1=h/h0-1, where h0=m_0*c^2,
 P   - pressure (in dyn/cm^2),''/
 Gamma - adiabatic exponent, Gamma= d ln P / d ln n.
 If MODE=1,2,or 3 then an input argument is XN,RHO, or H1, respectively;
    the other 4 arguments are fitted on output. */
{
   const double c2=8.9875518e20 /* c^2 in CGS */ ;
   double RLG,PLG,X,Y,Z; 
      if (MODE==1) {
         X=*XN;
         Y=FitRofN(KEOS,X);
         *RHO=Y;
      } else if (MODE==2) {
         X=*RHO;
         Y=FitNofR(KEOS,X);
         *XN=Y;
      } else if (MODE==3) {
         X=*H1;
         RLG=FitHEOS(KEOS,X);
         *RHO=pow(10.,RLG);
      } else {
         printf("NSEOSFIT: unknown MODE\n");
         exit(0);
      }
      if (KEOS<1 || KEOS>2) {
         printf("NSEOSFIT: unknown EOS\n");
         exit(0);
      }
      X=*RHO;
      RLG=log10(X);
      Y=FitGammaP(KEOS,RLG,&PLG);
      *Gamma=Y;
      Y=pow(10.,PLG)+3.5e14*X;
      *P=Y;
      if (MODE==3) {
         Z=*H1;
         Y=(X+Y/c2)/(Z+1.)/1.66e15;
         *XN=Y;
      } else {
         Z=*XN;
         Y=(X+Y/c2)/Z/1.66e15-1.;
         *H1=Y;
      }
   return 0;
}
/******************** END OF NSEOSFIT ******************/

double FitRofN(int KEOS, double XN)
/*                                                       Version 31.08.04
 Fit of (rho(n_b)-m_0*n_b) consistent with FitEOS to within
 2.1/5.4% for SLy4, 2.0/6.6% for FPS,
 Error in rho within 0.47% for SLy4, 0.35% for FPS, 0.8% for SLy+APR
 Input: KEOS; XN=n_b(fm^{-3})
 Output: RHO = rho[g/cc]
 KEOS=1: SLy, KEOS=2: FPS
*/
{
   static double A[3][7] = {{.423,2.42,.031,.78,.238,.912,3.674},
            {.32,2.17,.173,3.01,.54,.847,3.581},
            {.342,2.23,2.2e-6,10.92,0.,.839,3.68}};
   int K;
   double RHO, F1, F2, XLG, FRM, F;
      if (KEOS<1 || KEOS>2) {printf("FitEOS: KEOS\n"); exit(0);};
      if (XN<1.e-16) {printf("FitRofN: too low XN\n"); exit(0);};
      if (XN>5.) {printf("FitRofN: too high XN\n"); exit(0);};
      K=KEOS-1;
      F1=(A[K][0]*pow(XN,A[K][1])+A[K][2]*pow(XN,A[K][3]))/pow(1.e0+A[K][4]*XN,2);
      F2=XN/(8.e-6+pow(XN,.585)*2.1);
      XLG=log10(XN);
      FRM=FERMI(A[K][5]*(XLG+A[K][6]));
      F=FRM*F2+(1.-FRM)*F1;
      RHO=XN*1.66e15*(1.+F);
   return RHO;
}
/******************** END OF FitRofN ******************/

double FitNofR(int KEOS, double RHO)
/*                                                       Version 30.08.04
 Fit of (rho-m_0*n_b(rho)) consistent with FitEOS to within
 2.1/5.4% for SLy4, 2.1/6.6% for FPS
 Error in n_b within .56% for SLy, .74% for FPS, 1.4% for SLy+APR
 Input: KEOS; RHO = rho[g/cc]
 Output: XN=n_b(fm^{-3})
 KEOS=1: SLy4, KEOS=2: FPS
*/
{
   static double A[3][7] = {{.183,1.26,6.88,3.612,2.248,.911,11.56},
            {.608,2.41,2.39,3.581,1.681,.850,11.64},
            {.173,1.18,9.97,3.787,2.634,.917,11.56}};
   int K;
   double X, F1, F2, RLG, FRM, F, XN;
      if (KEOS<1 || KEOS>2) {printf("FitEOSt: KEOS\n"); exit(0);};
      if (RHO<.1) {printf("FitNofR: too low RHO\n"); exit(0);};
      if (RHO>1.e16) {printf("FitNofR: too high RHO\n"); exit(0);};
      X=RHO/1.66e15;
      K=KEOS-1;
      F1=(A[K][0]*pow(X,A[K][1])+A[K][2]*pow(X,A[K][3]))/pow(1.+A[K][4]*X,3);
      F2=X/(8.e-6+2.1*pow(X,.585));
      RLG=log10(RHO);
      FRM=FERMI(A[K][5]*(RLG-A[K][6]));
      F=FRM*F2+(1.-FRM)*F1;
      XN=X/(1.+F);
   return XN;
}
/******************** END OF FitNofR ******************/

double FitHEOS(int KEOS, double H1)
/*                                                       Version 30.08.04
 Fit of NS EOS for rotating configurations
 Input: KEOS; H1=H/H0-1, H0=H(RHO=7.86 g/cc) \approx 1
 Output: RLG = lg(rho[g/cc])
 KEOS=1: SLy4, rho=1.e4-4.e15 g/cc,
    average error 1.3%, max.error 4.0% at H1=0.0095 (rho=3.69e11 g/cc)
 KEOS=2: FPS, rho=1.e4-1.e17 g/cc,
    average error 0.95%, max.error 4.2% at H1=0.0225 (rho=1.23e14 g/cc)
*/
{
   static double A[3][16] = {{5.926,.4704,20.13,.2347,3.07,97.8,-2.012,89.85,
        34.96,15.328,.621,63.1,68.5,2.518,2.6,1.363},
            {5.926,.4704,19.92,.2333,2.63,54.7,-1.926,36.89,
       11.97,15.432,.6731,49.40,11.47,1.425,3.0,.913},
            {5.926,.4704,19.927,.2332,2.90,90.2,-2.002,37.254,
       11.95,15.259,.5080,31.18,11.35,.678,3.0,.762}};
   int K;
   double X, RLG1, RLG3, G, RLG2, C1, C2, RLG;
      if (KEOS<1 && KEOS>2) {printf("FitEOS: KEOS\n"); exit(0);};
      if (H1<0.) {printf("FitHEOS: negative H1\n"); exit(0);};
      if (H1>5.) {printf("FitHEOS: too high H1\n"); exit(0);};
      X=log10(H1);
      K=KEOS-1;
      RLG1=A[K][0]+A[K][1]*X+A[K][2]*pow(H1,A[K][3])/(1.+A[K][4]*H1);
      RLG3=(A[K][9]+A[K][10]*X);
      G=pow(A[K][11]*H1,7);
      RLG2=(A[K][7]+A[K][8]*X+G*RLG3)/(1.+A[K][12]*H1+G);
      C1=A[K][5];
      C2=A[K][6];
      RLG=RLG1*FERMI(C1*(X-C2))+RLG2*FERMI(C1*(C2-X))+
       FERMI(A[K][14]*(A[K][15]-X))*A[K][13];
  return RLG;
}
/******************** END OF FitHEOS ******************/

double FERMI(double X)
{
   double F;
      if (X>40.) {return 0.;};
      if (X<-40.) {return 1.;};
      F=1./(exp(X)+1.);
   return F;
}
/******************** END OF FERMI ******************/

double FitGammaP(int KEOS, double RLG, double *PLG)
/*                                                       Version 30.08.04
 Adiabatic index - formula consistent with FitEOS
 Input: KEOS; RLG = lg(rho[g/cc])
 Output: PLG=lg(P[dyn/cm^2]), Gamma=d log P / d log n_b
 KEOS=1: SLy4
 KEOS=2: FPS
*/
{
  const double c2=pow(2.99792458e10,2);
  static double A[3][18] = {{6.22,6.121,.005925,.16326,6.48,11.4971,19.105,.8938,
     6.54,11.4950,-22.775,1.5707,4.3,14.08,27.80,-1.653,1.50,14.67},
            {6.22,6.121,.006004,.16345,6.50,11.8440,17.24,1.065,
     6.54,11.8421,-22.003,1.5552,9.3,14.19,23.73,-1.508,1.79,15.13},
            {6.22,6.121,.006035,.16354,4.73,11.5831,12.589,1.4365,
     4.75,11.5756,-42.489,3.8175,2.3,14.81,29.80,-2.976,1.99,14.93}};
  int K;
  double X, Y, D4, P01, P02, P03, P04, F1,F2, F3, F4, DZDX, Gamma;
      if (KEOS<1 || KEOS>2) {printf("FitEOS: KEOS\n"); exit(0);};
      if (RLG<-1.) {printf("FitGammaP: too low RLG\n"); exit(0);};
      if (RLG>16.) {printf("FitGammaP: too high RLG\n"); exit(0);};
      X=RLG;
      K=KEOS-1;
      D4=(1.+A[K][3]*X);
      P01=(A[K][0]+A[K][1]*X+A[K][2]*pow(X,3))/D4;
      P02=A[K][6]+A[K][7]*X;
      P03=A[K][10]+A[K][11]*X;
      P04=A[K][14]+A[K][15]*X;
      F1=FERMI(A[K][4]*(X-A[K][5]));
      F2=FERMI(A[K][8]*(A[K][9]-X));
      F3=FERMI(A[K][12]*(A[K][13]-X));
      F4=FERMI(A[K][16]*(A[K][17]-X));
      Y=P01*F1+P02*F2+P03*F3+P04*F4;
      *PLG=Y;
      DZDX=F1/D4               /* d\zeta/dX=d ln P/d ln rho */
     *  ((A[K][1]-A[K][0]*A[K][3]+3.*A[K][2]*X*X
     +    2.*A[K][2]*A[K][3]*pow(X,3))/D4
     -  A[K][4]*(A[K][0]+A[K][1]*X+A[K][2]*pow(X,3))
     *    FERMI(A[K][4]*(A[K][5]-X)))
     +  F2*(A[K][7]+A[K][8]*FERMI(A[K][8]*(X-A[K][9]))*P02)
     +F3*(A[K][11]+A[K][12]*FERMI(A[K][12]*(X-A[K][13]))*P03)
     +F4*(A[K][15]+A[K][16]*FERMI(A[K][16]*(X-A[K][17]))*P04);
      Gamma=(1.+pow(10.,Y-RLG)/c2)*DZDX;
      return Gamma;
}
/******************** END OF FitGammaP ******************/