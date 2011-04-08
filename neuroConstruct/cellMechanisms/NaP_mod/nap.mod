TITLE nap-channel
COMMENT
	Based on D'Angelo (2001)
ENDCOMMENT
 
NEURON { 
	SUFFIX %Name%
	USEION na READ ena WRITE ina 
	RANGE gmax, ina
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
    gmax= %Max Conductance Density% (mho/cm2): (PG: was 0.00002) after rescaling, at T=23


	Aalpha_m = 0.091 (/ms-mV):in case of sig/expn this is /ms
	Kalpha_m = 5(mV) : 
	V0alpha_m = -42(mV)
	
	Abeta_m = -0.062(/ms)
	Kbeta_m = -5(mV): R
	V0beta_m = -42(mV)

		
	ena = 87.4(mV) 
	celsius = 30 (degC) 
	Q10 = 3 :
	TC =1
} 

STATE { 
	m 
	
} 

ASSIGNED { 
	v (mV) 
	ina (mA/cm2) 
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
	ina = g*(v - ena)
	alpha_m = alp_m(v)
	beta_m = bet_m(v) 

    : PG for testing/// printf("orig mod -- t=%g v=%g m=%g alpha_m=%g beta_m=%g \n ", t, v, m, alpha_m, beta_m)
} 
 
DERIVATIVE states { 
	rate(v) 
	m' =(m_inf - m)/tau_m 
	
} 
 
FUNCTION alp_m(v(mV))(/ms) {            
	alp_m = TC*Aalpha_m*linoid(v,V0alpha_m,Kalpha_m) 
} 
 
FUNCTION bet_m(v(mV))(/ms) {            
	bet_m = TC*Abeta_m*linoid(v,V0beta_m,Kbeta_m) 
} 
 
 
PROCEDURE rate(v (mV)) {LOCAL a_m, b_m
	TABLE m_inf, tau_m
	DEPEND Aalpha_m, Kalpha_m, V0alpha_m, 
	       Abeta_m, Kbeta_m, V0beta_m,
               celsius FROM -100 TO 100 WITH 200 
	a_m = alp_m(v)  
	b_m = bet_m(v) 
	m_inf = sig(v,-42,5): Actually -42 (mV),etc..
	tau_m = 5/(a_m + b_m) 
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

: Special infinity values still to be implemented
