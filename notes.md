# Notes

## TODO

### Pfxt, Sfxt, test, dump, cppr

## Notations

### res = resistance

### PfxtNode = ?, Sfxt = ?

### el = early late

### rf = rise fall

### frf = from pin rise fall

### trf = to pin rise fall

### at = arrival time

### rat = require time

### _fprop_at = forward propagation arrival time

### _bprop_rat = backward propagation require time

#### \*\_delay\[el][frf][trf]: _delay is a lut,\*_delay\[el][frf][trf] returns a single float value after looking-up

#### _rat\[fel][frf].emplace(arc, tel, trf, val)

***

## System Requirements

### It is developed completely from the ground up using C++17

### GNU C++ Compiler v7.3 with -std=c++1z

### tclsh

### cmake

***

## Licenses

### The core of OpenTimer is under MIT license. However, it is built, tested, and documented using several third-party tools and services. Thanks a lot

### units: a compile-time header-only dimensional analysis and unit conversion

### Parser-SPEF: a fast C++ header-only parser for standard parasitic exchange format (SPEF)

### PEGTL: parsing expression grammar template library

### Cpp-Taskflow: fast C++ parallel programming with task dependencies

### NanGate 45nm Library: open-source standard-cell library for testing and exploring EDA flows

### OSU PDK: Oklahoma State University system on chip (SoC) design flows

### Synopsys TAP-in: Synopsys technology access program for liberty user guide and open-source SDC parser
