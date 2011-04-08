TITLE ksl-channel
COMMENT
	Based on D'Angelo(2001)
ENDCOMMENT
 
NEURON { 
	SUFFIX %Name%
	USEION k READ ek WRITE ik 
	RANGE gmax, ik
	RANGE g, alpha_m, beta_m
	RANGE m_inf, tau_m
	
	RANGE Aalpha_m, Kalpha_m, V0alpha_m
	RANGE Abeta_m, Kbeta_m, V0beta_m
	GLOBAL TC, Q10
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 

    gmax= %Max Conductance Density% (mho/cm2)  :   PG: was 0.00035, after rescaling, at T=23

	Aalpha_m = 0.008(/ms):in case of sig/expn this is /ms
	Kalpha_m = 40(mV) : 
	V0alpha_m = -30(mV)
	
	Abeta_m = 0.008(/ms)
	Kbeta_m = -20(mV): R
	V0beta_m = -30(mV)
	
	Km_inf = 6 (mV)
	V0m_inf = -30 (mV)
		
	ek = -84.7(mV) 
	celsius = 30 (degC) 
	Q10 = 3 :
	TC =1
} 

STATE { 
	m 
} 

ASSIGNED { 
	v (mV) 
	ik (mA/cm2) 
	m_inf 
	tau_m (ms) 
	g (mho/cm2) 
	alpha_m (/ms)
	beta_m (/ms)
} 
 
INITIAL { 
	rate(v) 
	m = m_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD cnexp 
	g = gmax*m
	ik = g*(v - ek)
	alpha_m = alp_m(v)
	beta_m = bet_m(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	m' =(m_inf - m)/tau_m 
} 
 
FUNCTION alp_m(v(mV))(/ms) {            
	alp_m = TC*Aalpha_m*expn(v,V0alpha_m,Kalpha_m) 
} 
 
FUNCTION bet_m(v(mV))(/ms) {            
	bet_m = TC*Abeta_m*expn(v,V0beta_m,Kbeta_m) 
} 
 
 
PROCEDURE rate(v (mV)) {LOCAL a_m, b_m
	TABLE m_inf, tau_m
	DEPEND Aalpha_m, Kalpha_m, V0alpha_m, 
	       Abeta_m, Kbeta_m, V0beta_m,
               celsius FROM -100 TO 100 WITH 200 
	a_m = alp_m(v)  
	b_m = bet_m(v) 
	
	m_inf = sig(v,V0m_inf,Km_inf)
	tau_m = 1/(a_m + b_m) 
} 

FUNCTION linoid(V (mV),x (mV),y (mV)) (mV) {
        if (fabs((V-x)/y) < 1e-6) {
                linoid = y*(1 + (V-x)/(2*y))
        } else {
	  linoid = (V-x)/(1-exp(-(V-x)/y))
	}
}

FUNCTION expn(V (mV), x (mV), y (mV)) {
	expn = exp((V-x)/y)
}

FUNCTION sig(V (mV),x (mV),y (mV)) {
	sig = 1/(1+exp(-(V-x)/y))
}

: Special infinity values have still to be implemented
