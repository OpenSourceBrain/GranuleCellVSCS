: kir_ro channel (potassium inward rectifier)

COMMENT
Based upon Table 2 of Rossi_JNS_98
But with increased gbar
ENDCOMMENT

NEURON {
  SUFFIX %Name%
  USEION k READ ek WRITE ik
  RANGE gmax, ik
  RANGE Aalpha_m, Kalpha_m, V0alpha_m
  RANGE Abeta_m, Kbeta_m, V0beta_m
  RANGE shift
  }

UNITS {
  (molar) = (1/liter)
  (mM)    = (millimolar)
  (mV)    = (millivolt)
  (mA)    = (milliamp)
  (S)     = (siemens)
  FARADAY = (faraday) (kilocoulombs)
  R       = (k-mole) (joule/degC)
  }

PARAMETER {
  gmax = %Max Conductance Density% (S/cm2): PG: was 9.5e-4

  shift = 0 (mV): DA 0

  Aalpha_m = 0.4 (/ms): 
  Kalpha_m = -0.08: ~24.4 :
  V0alpha_m = 86.3(mV) : 

  Abeta_m = 0.4 (/ms): 
  Kbeta_m = 0.04 (mV): 
  V0beta_m = 86.3(mV): 
  }

ASSIGNED {
  celsius (degC) : Simulation should be at 30
  v        (mV)
  ek       (mV)
  ik       (mA/cm2)  : or should it be specified as the f(ast) part?
  }

STATE {
  m 
  }

BREAKPOINT {
  SOLVE states METHOD cnexp
  ik = gmax * m*(v-ek)

    :printf("orig mod -- t=%g v=%g m=%g\n", t, v, m)
  }

DERIVATIVE states {
  m' = malpha(v)*(1-m)-mbeta(v)*m
  }

INITIAL {
  m = malpha(v)/(malpha(v)+mbeta(v))
  }

FUNCTION malpha(v (mV)) (/ms) {
  UNITSOFF  
  malpha = Aalpha_m*exp(Kalpha_m*(v+V0alpha_m))
  UNITSON

  }

FUNCTION mbeta(v (mV)) (/ms) {
  UNITSOFF
  mbeta = Abeta_m*exp(Kbeta_m*(v+V0beta_m))
  UNITSON
  }
  
