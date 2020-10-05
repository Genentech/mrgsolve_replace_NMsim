;; 1. Based on: 
;; 2. Description: Model 6 - 2cmt model with covariates
;; x1. Author: user

; Reference: 
; Nonlinear Mixed‐Effects Model Development and Simulation Using nlmixr and Related R Open‐Source Packages
; https://doi.org/10.1002/psp4.12445

$PROB examplomycin

$INPUT ID TIME DV CMT WT SEX AMT EVID PHASE=DROP

$DATA ../data/examplomycin2.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9 

$MODEL NCOMP=3

$PK
NMID = ID

; Covariates
WT2 = WT
IF(WT.LE.0) WT2 = 70


MU_1 = LOG(THETA(1)) + THETA(11) * LOG(WT2/70)
CL   = EXP(MU_1 + ETA(1))

MU_2 = LOG(THETA(2)) + THETA(12) * SEX
V2   = EXP(MU_2 + ETA(2))

MU_3 = LOG(THETA(3)) 
V3   = EXP(MU_3 + ETA(3))

MU_4 = LOG(THETA(4)) 
Q    = EXP(MU_4 + ETA(4))

MU_5 = LOG(THETA(5)) 
KA   = EXP(MU_5 + ETA(5))

K   = CL/V2
K23 = Q /V2	
K32 = Q /V3

; Calc TAD
;; initialize for new individual or reset record
IF (NEWIND.LT.2.OR.EVID.EQ.3) THEN
 TDOS=-9999 ; Time of most recent dose. -9999 if no previous dose.
 TAD=0 ; Time After Dose
ENDIF
IF (EVID.EQ.1.OR.EVID.EQ.4) TDOS=TIME
IF (TDOS.GT.-9999) TAD=TIME-TDOS


$DES
DADT(1) =-KA*A(1)                           ; Depot
DADT(2) = KA*A(1) - (K+K23)*A(2) + K32*A(3) ; Central
DADT(3) =              K23 *A(2) - K32*A(3) ; Peripheral


$ERROR
IPRED = (A(2)/V2+.000001)

PROP = IPRED*THETA(6)
ADD  = THETA(7)
SD   = SQRT(PROP*PROP+ADD*ADD)

Y = IPRED + EPS(1)*SD
IRES = DV - IPRED


$THETA   
(0, 0.1)  ; 1 CL
(0, 2)    ; 2 V2
(0, 5)    ; 3 V3
(0, 0.3)  ; 4 Q
(0, 1)    ; 5 KA
(0, 0.1)  ; 6 Prop err
0 FIX     ; 7 Add err
0 FIX  ; 8 blank
0 FIX  ; 9 blank
0 FIX  ; 10 blank
(0.75)  ; 11 WT on CL
(-0.2) ; 12 SEX on V2

$OMEGA BLOCK(2)
.1 ; 1 CL 
.03 .1 ; 2 V2 
$OMEGA 
.1 ; 3 V3 
.1 ; 4 Q 
.1 ; 5 KA
$SIGMA 1 FIX

$ESTIMATION METHOD=COND INTERACTION MAXEVAL=9999 NSIG=3 SIGL=10
PRINT=1 NOTHETABOUNDTEST NOABORT
$COVARIANCE PRINT=E UNCONDITIONAL SIGL=10

$TABLE ID NMID TIME TAD DV CMT WT SEX AMT EVID ETAS(1:LAST) CL V2
IPRED CWRES IRES NOPRINT ONEHEADER FILE=run006.tab

