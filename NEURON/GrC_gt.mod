TITLE GrC tonic conductance

UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
	(uS) = (microsiemens)
} 

NEURON {
	POINT_PROCESS GrC_gt
	NONSPECIFIC_CURRENT i
	RANGE gmax, e, i, stdel, stdur
}

PARAMETER {
	gmax = 0.0 (uS)
	e = 0.0 (mV)
	stdel = 0.0 (ms)
	stdur = 0.0 (ms)
}

ASSIGNED {
	v (mV)
	i (nA)
}

INITIAL {
	i = 0.0
}

BREAKPOINT { LOCAL stend
	stend = stdel + stdur
	if ((t >= stdel) && (t <= stend)) {
		i = gmax*(v - e)
	} else {
		i = 0.0
	}
}

