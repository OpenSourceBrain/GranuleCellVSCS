TITLE : cacm_ms10 (Changes to intracellular calcium concentration)

COMMENT
	Based on MS_model(98)
ENDCOMMENT

NEURON {
	SUFFIX %Name%
	USEION ca READ ica WRITE cai 
	RANGE depth, tau, cai0
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	 F   = (faraday) (coulombs)
}

PARAMETER {
	depth = 84     (nm)  : assume volume = area*depth ; depth=d(article)
	                     : assume nm is recognized
    tau = 10       (ms)  : assuming beta in article is (/ms)
	cai0 = 0.75e-4 (mM)  : Resting intracellular concentration ://it was 0.755
                         : Requires explicit use in INITIAL
			             : block for it to take precedence over cai0_ca_ion
			             : Do not forget to initialize in hoc if different
			             : from this default.
      }

ASSIGNED {
	ica (mA/cm2)
}

STATE {
	cai (mM)
}

INITIAL {
	cai = cai0
}

BREAKPOINT {
	SOLVE integrate METHOD derivimplicit : Why derivimplicit?
}

DERIVATIVE integrate {
	cai' = -(1.0e7)*ica/(2*F*depth) - (cai - cai0)/tau
}
