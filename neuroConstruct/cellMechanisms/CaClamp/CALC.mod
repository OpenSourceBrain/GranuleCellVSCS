
TITLE Calcium Clamp

COMMENT

        Calcium "Clamp": not really physiological
        Used to clamp the Ca2+ concentration to a higher level for a time period.
        Useful for testing Ca2+ dependent channels
        
ENDCOMMENT

NEURON {
        SUFFIX %Name%
        USEION ca READ ica WRITE cai
        RANGE del, dur, level1, level2, targetLevel, rate
}

UNITS {
        (mV)    = (millivolt)
        (mA)    = (milliamp)
	(um)    = (micron)
	(molar) = (1/liter)
        (mM)    = (millimolar)
}

PARAMETER {
        ica             (mA/cm2)
        del = 100           (ms)
        dur = 20           (ms)
        level1 = 5e-5     (mM)
        level2 = 5e-4     (mM)
        rate = 1
}

ASSIGNED {
	ca_pump_i	(mA)
	targetLevel	(mM)
}
STATE {
	cai (mM)
}

INITIAL {
        cai = level1
        targetLevel = level1
}

BREAKPOINT {

        SOLVE conc METHOD derivimplicit
}

DERIVATIVE conc {
        if (t > del && t<= (del+dur)) {
            targetLevel = level2
        }
        else {
            targetLevel = level1
        }

	cai' =  rate * (targetLevel - cai)
}


