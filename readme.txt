========================================================================
    asteroid Project Overview
========================================================================

// Author Stephen Bryan
// June 22, 2017

Project to simulate gravitational forces on N bodies (data contained in initialization file), 
and expecially collisions among those bodies.

Initially, model collisions as elastic, but make more realistic down the road.

Usage:  asteroid.exe [filename]

If filename is not specified, "init.data" is looked for

Format of data file:

Options: [momtol = 2e-2,] [dvnormtol = 2e-2,] [runs = 1.0e9,] [tstep = 1.0,] [minvforcheck = 0.01,] [startdatasave=2e6,] [nthdatasave=1024]
<Object data>*

All options are optional.  General format is 
    optnam = val,
with a final comma unnecessary.  Whitespace is ignored, the equals sign and commas as separators are required.

Current options and defaults are:
  momtol = 0.01 (max momentum change before timestep will be shortened)
  dvnormtol = 0.01 (max momentum change in the normal diection before timestep will be shortened)
  runs = 1.0e6 (number of nominal timesteps)
  tstep = 1.0 (nominal timestep)
  minvforcheck = 0.01 (minimum velocity below which momentum constraint momtol will not be checked)
  startdatasave=0 (specify the time to start saving data - used to ignore behavior before collisions)
  nthdatasave=1 (save only every Nth data point)

<Object data>
xloc yloc zloc xvel yvel zvel mass

Any number of objects can be specified, each to its own line.


Project TBD:
- write, calculate, save validation tests
- add randomness to rho in Object::Object, or in Object data (initialization file)
- enhance visualization
- generate object data with user-specified characteristics
- test 
	- simple tests for validation
	- larger number of objects around massive center
	- binary pair at center
	- objects with no massive center
- make collisions more realistic
- implement accretion






/////////////////////////////////////////////////////////////////////////////
