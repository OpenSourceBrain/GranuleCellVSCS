TITLE kv-channel kv_ms10.mod
COMMENT
	Based on Maex/De Schutter
ENDCOMMENT
 
NEURON { 
	SUFFIX %Name%
	USEION k READ ek WRITE ik
	RANGE gmax, ik
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

    gmax= %Max Conductance Density% (mho/cm2): (PG: was 0.003) after rescaling, at T=22

	shift = -10(mV)

	Aalpha_m = 0.17(/ms):in case of sig/expn this is /ms
	Kalpha_m = 13.7(mV) :
	V0alpha_m = -38(mV)
	
	Abeta_m = 0.17(/ms)
	Kbeta_m = -55.55(mV): 
	V0beta_m = -38(mV)

	Aalpha_h = 6.0e-5(/ms)
	Kalpha_h = -12.5(mV) : R=
	V0alpha_h = -46 (mV)
 
	Abeta_h = 1.1e-3(/ms)
	Kbeta_h = 12.4(mV)
	V0beta_h = -44(mV) : or should it be +44?
	
	ek = -90 (mV) 
	celsius = 37 (degC) 
	Q10 = 3 :
	TC=5
} 

STATE { 
	m 
	h 
} 

ASSIGNED { 
	v (mV) 
	ik (mA/cm2) 
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
	g = gmax*m*m*m*m*h 
	ik = g*(v - ek)
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
	alp_h = TC*(0.0007 (/ms) + Aalpha_h*expn(v+shift,V0alpha_h,Kalpha_h))
	if (alp_h < TC*7.6e-4){                                              :::::::::::::::::::::::::::   corrected by PG, origially was if (alp_h < 7.6e-4){
	alp_h = TC*7.6e-4
	} 
} 
 
FUNCTION bet_h(v(mV))(/ms) {             
	bet_h = TC*Abeta_h*sig(v+shift,V0beta_h,Kbeta_h)
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


