TITLE kca-channel kca_ms10.mod
COMMENT
	Based on MS(98)
ENDCOMMENT
 
NEURON { 
	SUFFIX %Name%
	USEION k READ ek WRITE ik
	USEION ca READ cai VALENCE 2
	RANGE gmax, ik
	RANGE g, alpha_m, beta_m
	RANGE m_inf, tau_m
	
	RANGE Aalpha_m, Kalpha_m, V0alpha_m
	RANGE Abeta_m, Kbeta_m, V0beta_m
	GLOBAL TC,Q10, shift
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	cai      (mM)  : or in ASSIGNED block? In TN also here
	shift  = -10(mV)
	
	gmax= %Max Conductance Density% (mho/cm2): (PG: this was 0.018) after rescaling, at T=22
	ek = -90 (mV) 
	celsius = 37 (degC) 
	Q10 = 3 :
	TC=5
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
	UNITSOFF  
 	alp_m = TC*2.5/(1+1.5e-3 * exp(-0.085*(v+shift))/cai) 
  	UNITSON                                    
} 
 
FUNCTION bet_m(v(mV))(/ms) {
 	UNITSOFF             
	bet_m = TC*1.5/(1+ cai/(1.5e-4 *exp(-0.077*(v+shift)))) 
      UNITSON
}
 
PROCEDURE rate(v (mV)) {LOCAL a_m, b_m
	a_m = alp_m(v)  
	b_m = bet_m(v) 
	m_inf = a_m/(a_m + b_m) 
	tau_m = 1/(a_m + b_m) 	
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
: Parameters should be made accessible from hoc
: Use of cac instead of ca to prevent any misunderstanding