: NetStim with flat distribution of output events
: spike number becomes meaningless, use duration instead

NEURON  {
  ARTIFICIAL_CELL NetStimF
  RANGE y
  RANGE interval, number, start
  RANGE noise, duration
}

PARAMETER {
        interval = 10 (ms) <1e-9,1e9>	: time between spikes (msec)
        number  = 10 <0,1e9>    	: number of spikes
        start = 50 (ms)       		: start of first spike
        noise = 0 <0,1>       		: amount of randomeaness (0.0 - 1.0)
	duration = 1000 (ms)		: burst duration
}

ASSIGNED {
        y
        event (ms)
        on
        end (ms)
}

PROCEDURE seed(x) {
        set_seed(x)
}

INITIAL {
        on = 0
        y = 0
        if (noise < 0) {
                noise = 0
        }
        if (noise > 1) {
                noise = 1
        }
        if (start >= 0 && number > 0) {
                : randomize the first spike so on average it occurs at
                : start + noise*interval
                event = start + invl(interval) - interval*(1. - noise)
                : but not earlier than 0
                if (event < 0) {
                        event = 0
                }
                net_send(event, 3)
        }
}

PROCEDURE init_sequence(t(ms)) {
        if (number > 0) {
                on = 1
                event = t
                end = t + 1e-6 + duration
        }
}

FUNCTION invl(mean (ms)) (ms) {
        if (mean <= 0.) {
                mean = .01 (ms) : I would worry if it were 0.
        }
        if (noise == 0) {
                invl = mean
        }else{
                invl = (1. - noise)*mean + noise*mean*exprand(1)
        }
}

PROCEDURE event_time() {
        if (number > 0) {
                event = event + invl(interval)
        }
        if (event > end) {
                on = 0
        }
}

NET_RECEIVE (w) {
        if (flag == 0) { : external event
                if (w > 0 && on == 0) { : turn on spike sequence
                        init_sequence(t)
                        net_send(0, 1)
                }else if (w < 0 && on == 1) { : turn off spiking
                        on = 0
                }
        }
        if (flag == 3) { : from INITIAL
                if (on == 0) {
                        init_sequence(t)
                        net_send(0, 1)
                }
        }
        if (flag == 1 && on == 1) {
                y = 2
                net_event(t)
                event_time()
                if (on == 1) {
                        net_send(event - t, 1)
                }
                net_send(.1, 2)
        }
        if (flag == 2) {
                y = 0
        }
}

