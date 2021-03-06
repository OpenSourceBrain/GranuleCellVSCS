TITLE ka-current ka_ms10.mod
COMMENT
	Based on MS(98)
ENDCOMMENT
 
NEURON { 
	SUFFIX ka_ms10
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	RANGE g
	RANGE mf, hf, tm, th
	
	RANGE Aalpha_m, Kalpha_m, V0alpha_m
	RANGE Abeta_m, Kbeta_m, V0beta_m
	
	RANGE Abeta_h, Kbeta_h, V0beta_h
	GLOBAL TC,Q10, shift
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 

	:
	shift = -10 (mV)
	
	Aalpha_m = 0.082(/ms):in case of sig/expn this is /ms
	Kalpha_m = -42.8(mV) :
	V0alpha_m = -43.5(mV)
	
	Abeta_m = 0.2(/ms)
	Kbeta_m = 19.8(mV): 
	V0beta_m = -46.7(mV)
	
	Aalpha_h = 0 : dummies (they are actually not used)
	Kalpha_h = 0
	V0alpha_h = 0

	
 
	Abeta_h = 0.2(/ms)
	Kbeta_h = -8.4(mV)
	V0beta_h = -78.8(mV) 
	
	
	gkbar= 0.0011  (mho/cm2): after rescaling, at T=23  
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
	
	mf
	hf
	tm (ms)
	th (ms)
	g (mho/cm2) 
} 
 
INITIAL { 
	rate(v) 
	m = mf
	h = hf 
	
} 
 
BREAKPOINT { 
	SOLVE states METHOD cnexp
 
	g = gkbar*m*m*m*h 
	ik = g*(v - ek)
} : extra variables for readout 
 
DERIVATIVE states { 
	rate(v) 
	m' =(mf- m)/tm 
	h' =(hf- h)/th
} 
 
FUNCTION tau_m(v(mV))(/ms) {            
	tau_m = TC*(Aalpha_m*expn(v+shift,V0alpha_m,Kalpha_m)+0.0334 (ms)) 
} 
 
FUNCTION m_inf(v(mV))(/ms) {            
	m_inf= TC*Abeta_m*sig(v+shift,V0beta_m,Kbeta_m) 
} 
 
FUNCTION tau_h(v(mV))(/ms) {             
	UNITSOFF	
	tau_h = TC*(2.16+0.006*v+(0.2/(57.9*exp(0.127*v)+0.000134*exp(-0.059*v))))
	UNITSON
} 
 
FUNCTION h_inf(v(mV))(/ms) {             
	h_inf= TC*Abeta_h*sig(v+shift,V0beta_h,Kbeta_h)
} 
 
PROCEDURE rate(v (mV)) {
	TABLE mf, hf, tm, th
	DEPEND Aalpha_m, Kalpha_m, V0alpha_m, 
	       Abeta_m, Kbeta_m, V0beta_m,
               
               Abeta_h, Kbeta_h, V0beta_h, celsius FROM -100 TO 100 WITH 200 
	tm= tau_m(v)  
	mf = m_inf(v) 
	th = tau_h(v)  
	hf = h_inf(v) 
	
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


