[PROB]
Source: mrgsolve model library

Two-compartment model for oral dosing 

Comment:  
Time unit: hour
Volume units: L

[GLOBAL]
#define C2 (CENT / V2) 

[CMT] @annotated
DEPOT  : Extravascular compartment
CENT   : Central compartment
PERIPH : Peripheral compartment (mass) 
AUCCENT: AUC in the central compartment

[PARAM] @annotated
TVCL   :  0.1 : Clearance (L/hour)
TVV2   :  2    : Central volume of distribution (L)
TVV3   :  5    : Peripheral volume of distribution (L)
TVQ    :  0.3  : Inter-compartmental clearnace (L/hour)
TVKA   :  1    : Absorption rate constant (1/hour)
TVALAG1:  0    : Lag time - EV dose
WT     : 70    : Individual weight (kg)
WTref  : 70    : Reference body weight (kg)
WTCL   :  0.75 : Effect of WT on CL
SEX    :  0    : 1 = male, 0 = female
SEXV2  : -0.2  : Effect of SEX=1 (male) on central compartment volume
PROP   :  0.1  : Proportional error (SD)
ADD    :  0    : Additive error (SD)

// ETA for sim with EBE
E_CLi : 0 : ETA on CL 
E_V2i : 0 : ETA on V2
E_V3i : 0 : ETA on V3
E_Qi  : 0 : ETA on Q
E_KAi : 0 : ETA on KA
  
[OMEGA] @annotated @block @name group1
// note: if you want to enter off-diagonal elements as correlations, add "@correlation" to this block
E_CL : 0.1       : ETA on CL
E_V2 : 0.03  0.1 : ETA on V2

[OMEGA] @annotated @name group2
E_V3   : 0.1  : ETA on V3
E_Q    : 0.1  : ETA on Q
E_KA   : 0.1  : ETA on KA

[SIGMA] 1

[MAIN]
// covariate effects on parameters
double CLCOV = pow(WT/WTref, WTCL);
double V2COV = exp(SEXV2 * SEX);

double CL   = TVCL * CLCOV * exp(E_CL + E_CLi);
double V2   = TVV2 * V2COV * exp(E_V2 + E_V2i); 
double V3   = TVV3 * exp(E_V3 + E_V3i);
double Q    = TVQ  * exp(E_Q  + E_Qi);
double KA   = TVKA * exp(E_KA + E_KAi); 

double K   = CL/V2;
double K23 = Q/V2;
double K32 = Q/V3;

ALAG_DEPOT = TVALAG1;

[ODE]
dxdt_DEPOT  =-KA*DEPOT;
dxdt_CENT   = KA*DEPOT - (K+K23)*CENT + K32*PERIPH;
dxdt_PERIPH =               K23 *CENT - K32*PERIPH;

dxdt_AUCCENT= C2;

[TABLE]

double IPRED = C2;
double SD = sqrt(PROP*PROP*IPRED*IPRED + ADD*ADD);
double DV = IPRED + EPS(1) * SD;

[CAPTURE] @annotated
IPRED : Concentration without residual variability
DV    : Concentration with residual variability 
  
[CAPTURE] E_CL E_V2 E_V3 E_Q E_KA CL V2

