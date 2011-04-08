TITLE ca-channel ca_ms10.mod
COMMENT
	Based on MS(98)
ENDCOMMENT
 
NEURON { 
	SUFFIX %Name%
	USEION ca READ eca WRITE ica
	RANGE gmax, ica
	RANGE g, alpha_m, beta_m, alpha_h, beta_h 
	RANGE m_inf, tau_m, h_inf, tau_h 
	
	RANGE Aalpha_m, Kalpha_m, V0alpha_m
	RANGE Abeta_m, Kbeta_m, V0beta_m
	RANGE Aalpha_h, Kalpha_h, V0alpha_h
	RANGE Abeta_h, Kbeta_h, V0beta_h
	GLOBAL TC,Q10,shift
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 

    gmax= %Max Conductance Density% (mho/cm2): (PG: was 0.0009) after rescaling, at T=23

	shift  = -10(mV)
	
	Aalpha_m = 1.6(/ms)
	Kalpha_m = 13.9(mV) 
	V0alpha_m = 5(mV)
	
	Abeta_m = -0.02(/ms-mV)
	Kbeta_m = -5(mV)
	V0beta_m = -8.9(mV)

	Aalpha_h = 0.005(/ms)
	Kalpha_h = -20(mV) 
	V0alpha_h = -60 (mV)
 
	Abeta_h = -0.005(/ms)
	Kbeta_h = -20(mV)
	V0beta_h = -60(mV) 
	
	eca = 80(mV) 
	celsius = 37 (degC) 
	Q10 = 3 
	TC=5
} 

STATE { 
	m 
	h 
} 

ASSIGNED { 
	v (mV) 
	ica (mA/cm2) 
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
	g = gmax*m*m*h 
	ica = g*(v - eca)
	alpha_m = alp_m(v)
	beta_m = bet_m(v) 
	alpha_h = alp_h(v)
	beta_h = bet_h(v) 

    :printf("-- mod, t: %g, v: %g, alpha_h: %g, beta: %g, tau: %g \n", t, v, alpha_h, beta_h, 1/(alpha_h + beta_h))
} 
 
DERIVATIVE states { 
	rate(v) 
	m' =(m_inf - m)/tau_m 
	h' =(h_inf - h)/tau_h 
} 
 
FUNCTION alp_m(v(mV))(/ms) {            
	alp_m = TC*Aalpha_m*sig(v+shift,V0alpha_m,Kalpha_m) 
} 
 
FUNCTION bet_m(v(mV))(/ms) {            
	bet_m = TC*Abeta_m*linoid(v+shift,V0beta_m,Kbeta_m) 
} 
 
FUNCTION alp_h(v(mV))(/ms) {             
	alp_h = TC*Aalpha_h*expn(v+shift,V0alpha_h,Kalpha_h) 
	if (alp_h > 0.005){
	alp_h = 0.005
	} : Consult RM on < or >
} 
 
FUNCTION bet_h(v(mV))(/ms) {             
	bet_h = TC*(0.005 (/ms)+Abeta_h*expn(v+shift,V0beta_h,Kbeta_h))
	if (bet_h < 0) {
	bet_h = 0
	}
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
	h_inf = a_h/(a_h + b_h) 
	tau_h = 1/(a_h + b_h) 
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
