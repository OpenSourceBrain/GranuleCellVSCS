TITLE naf-channel naf_ms_10.mod (shift=-10)
COMMENT
	MS_model:Based on Maex/DeSchutter (1998) and Gabbiani et al. (1994)
	Q: So V0alpha_m=-39 mV here (etc.), like in Gabbiani, so a shift towards -29
      mV still has to be made (denoted by shift = -10 mV)? Probably yes.
	Q: What means GLOBAL again, with regard to TC, Q10?
ENDCOMMENT
 
NEURON { 
	SUFFIX %Name%
	USEION na READ ena WRITE ina 
	RANGE gmax, ina
	RANGE g, alpha_m, beta_m, alpha_h, beta_h 
	RANGE m_inf, tau_m, h_inf, tau_h 
	
	RANGE Aalpha_m, Kalpha_m, V0alpha_m
	RANGE Abeta_m, Kbeta_m, V0beta_m
	RANGE Aalpha_h, Kalpha_h, V0alpha_h
	RANGE Abeta_h, Kbeta_h, V0beta_h
	GLOBAL TC, Q10, shift
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
    gmax= %Max Conductance Density%  (mho/cm2): (PG: was 0.0546) after rescaling, at T=22   

	shift = -10(mV)

	Aalpha_m = 1.5 (/ms):in case of sig/expn this is /ms
	Kalpha_m = 12.35(mV) : R=1/K=0.081
	V0alpha_m = -39(mV)
	
	Abeta_m = 1.5 (/ms)
	Kbeta_m = -15.15(mV): R=-0.066
	V0beta_m = -39(mV)

	Aalpha_h = 0.12 (/ms)
	Kalpha_h = -11.23(mV) : R=-0.089
	V0alpha_h = -50 (mV)
 
	Abeta_h = 0.12(/ms)
	Kbeta_h = 11.23(mV)
	V0beta_h = -50(mV)

	
	
	ena = 55(mV) 
	celsius = 37 (degC) 
	Q10 = 3 : not really used now (9-10-02)
	TC = 5
} 

STATE { 
	m 
	h 
} 

ASSIGNED { 
	v (mV) 
	ina (mA/cm2) 
	m_inf 
	h_inf 
	tau_m (ms) 
	tau_h (ms) 
	g (mho/cm2) 
	alpha_m (/ms)
	beta_m (/ms)
	alpha_h (/ms)
	beta_h (/ms)
      
} 
 
INITIAL { 
	rate(v) 
	m = m_inf 
	h = h_inf 
	
} 
 
BREAKPOINT { 
	SOLVE states METHOD cnexp  
	g = gmax*(m^3)*h 
	ina = g*(v - ena)
	alpha_m = alp_m(v)
	beta_m = bet_m(v) 
	alpha_h = alp_h(v)
	beta_h = bet_h(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	m' =(m_inf - m)/tau_m 
	h' =(h_inf - h)/tau_h 
} 
 
FUNCTION alp_m(v(mV))(/ms) {            
	alp_m = TC*Aalpha_m*expn(v+shift,V0alpha_m,Kalpha_m) 
} 
 
FUNCTION bet_m(v(mV))(/ms) {            
	bet_m = TC*Abeta_m*expn(v+shift,V0beta_m,Kbeta_m) 
} 
 
FUNCTION alp_h(v(mV))(/ms) {             
	alp_h = TC*Aalpha_h*expn(v+shift,V0alpha_h,Kalpha_h) 
} 
 
FUNCTION bet_h(v(mV))(/ms) {             
	bet_h = TC*Abeta_h*expn(v+shift,V0beta_h,Kbeta_h)
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_m, b_m, a_h, b_h 
	TABLE m_inf, tau_m, h_inf, tau_h 
	DEPEND Aalpha_m, Kalpha_m, V0alpha_m, 
	       Abeta_m, Kbeta_m, V0beta_m,
               Aalpha_h, Kalpha_h, V0alpha_h,
               Abeta_h, Kbeta_h, V0beta_h, celsius FROM -100 TO 100 WITH 200 
	a_m = alp_m(v)  
	b_m = bet_m(v) 
	a_h = alp_h(v)  
	b_h = bet_h(v) 
	m_inf = a_m/(a_m + b_m) 
	tau_m = 1/(a_m + b_m) 
	if (tau_m < 0.05/TC) {
	tau_m = 0.05/TC
	}
	h_inf = a_h/(a_h + b_h) 
	tau_h = 1/(a_h + b_h)
	if (tau_h < 0.225/TC) {
	tau_h = 0.225/TC
	} 
} 

FUNCTION linoid(V (mV),x (mV),y (mV)) (mV) {
        if (fabs((V-x)/y) < 1e-6) {
                linoid = y*(1 + (V-x)/(2*y))
        }else{
                linoid = (V-x)/(1-exp(-(V-x)/y))
        }
}

FUNCTION expn(V (mV),x (mV),y(mV)) {
	expn = exp((V-x)/y)
}

FUNCTION sig(V (mV),x (mV),y (mV)) {
	sig = 1/(1+exp(-(V-x)/y))
}

: The special Borg-Graham minimum time constants have to be implemented
