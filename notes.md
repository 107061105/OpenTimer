### el = early late
### rf = rise fall
### frf = from pin rise fall
### trf = to pin rise fall
### at = arrival time
### rat = require time
### _fprop_at = forward propagation arrival time
### _bprop_rat = backward propagation require time

#### *_delay[el][frf][trf]: _delay is a lut, *_delay[el][frf][trf] returns a single float value after looking-up
#### _rat[fel][frf].emplace(arc, tel, trf, val);