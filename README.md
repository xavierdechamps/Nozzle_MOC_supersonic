# Nozzle_MOC_supersonic
A code based on the Method of Characteristics to solve the supersonic flow inside a nozzle.

This Matlab / Octave code solves the flow inside a nozzle of known shape. The flow is assumed to be steady, two-dimensional, irrotational and supersonic. 
The methods developed in this project are strongly based on the reference book "Gas Dynamics, Volume II, Multidimensional Flow" by Zucrow Maurice J. and Hoffman Joe D., 
John Wiley and Sons, 1977.

The main routine is called **MOC_2D_steady_irrotational_main.m**. It contains all the inputs for the geometry of the nozzle. The geometry of this nozzle contains only 
the diverging part of the nozzle and a sonic condition is assumed at the throat. The radius of the throat is adapted by the parameter (geom.yt). The downstream region shows first 
a circular arc of radius (geom.rhod) which extends up to an angle (geom.ta). The diverging section is made of a parabolic curve which extends up to an axial distance (geom.xe) 
and shows an exit lip angle (geom.te). If (geom.ta) and (geom.te) are equal, then the diverging section is a line.

The initial-value line is chosen as the line where the y-velocity component is equal to zero. This line is discretized by (geom.NI) points and the data is propagated downstream 
by intersecting left-running and right-running characteristics. The compatibility and characteristic equations are numerically solved by a modified Euler predictor-corrector 
method to determine the x- and y-velocity components and the location of the new node.

![Example of intersections of characteristics](https://github.com/xavierdechamps/Nozzle_MOC_supersonic/blob/master/Images/Characteristics.jpg)

![Example of pressure distribution](https://github.com/xavierdechamps/Nozzle_MOC_supersonic/blob/master/Images/Pressure.jpg)

![Example of pressure distribution](https://github.com/xavierdechamps/Nozzle_MOC_supersonic/blob/master/Images/Mach_number.jpg)

