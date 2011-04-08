DEFINE nspikes 50

NEURON {
        POINT_PROCESS Bekkers3SynA
        NONSPECIFIC_CURRENT i
        RANGE taur, taud1, taud2, taud3, bamp, bpow, bfrac1, bfrac2, bfrac3
        RANGE i, g, gmax
}


UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (microsiemens)
}

PARAMETER {
	gmax = 0.001 (uS)
        bamp = 1.0
	bfrac1 = 2.23391 
	bfrac2 = 0.286137
	bfrac3 = 0.0830753 
	bpow = 11
	taur = 0.10337 (ms)
	taud1 = 0.451106 (ms)
	taud2 = 2.88445 (ms)
	taud3 = 21.6701 (ms)
        e = 0      (mV)
}


ASSIGNED {
        v (mV)
        i (nA)
        g (uS)
	tspike[nspikes] (ms)
}

INITIAL {
	LOCAL cspike
	cspike = 0
	while (cspike < nspikes) {
		tspike[cspike] = 0
		cspike = cspike + 1
	}
}


BREAKPOINT {
        g = gtrace(t)
        i = g*(v - e)
}


FUNCTION myexp(x) {
        if (x < -100) {
        	myexp = 0
        } else {
        	myexp = exp(x)
        }
}


FUNCTION gtrace(x) {
	LOCAL cspike
	cspike = 0
	gtrace = 0
	while ((cspike < nspikes) && (tspike[cspike] != 0)) {
       		gtrace = gtrace + (1 - myexp(-(x - tspike[cspike])/taur))^bpow * (bfrac1*myexp(-(x - tspike[cspike])/taud1) + bfrac2*myexp(-(x - tspike[cspike])/taud2) + bfrac3*myexp(-(x - tspike[cspike])/taud3))
		cspike = cspike + 1
	} 
	gtrace = gtrace*bamp*gmax
}


NET_RECEIVE(weight) {
	LOCAL cspike1
	cspike1 = nspikes - 1
	while (cspike1 > 0) {
        	tspike[cspike1] = tspike[(cspike1 - 1)]
		cspike1 = cspike1 - 1
	}
	tspike[0] = t
}
