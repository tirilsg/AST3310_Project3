# AST3310_Project3


# AST3310_Project2

This repository contains python code constructing a model for transportation of energy by convection in a 2D space. The model gives a more in depth look at both how fluids spread in a 2D space with boundaries, and how convectivley transported energy in the photosphere of the sun behaves. A report `AST_3310_Project3.pdf` can be assembled from the files `AST_3310_Project3.tex`, `references.bib` and the map `figures`, where the physics behind the model presented in `project3.py` is detailed, and results for relevant simulations presented and discussed. 

-------------------------------------

### `project3.py`: 
Contains a class `Convection2D`, which defnies a xy-space of size $12\cdot10^6\times4\cdot10^6$ [m]. The class initializes a temperature and pressure distribution in a hydrostatic equilibrium, and has a method which creates a temperature perturbation in the initial temperature condition which can be used to simulate behaviours. The class has methods; 

* `initialise():` which initialises the system
* `boundary_conditions():` which refedines the behaviour occurring at the boundaries of the 2D space
* `central_x(var):` a method for calculating derivatives by central differencing in the x-direction
* `central_y(var):` a method for calculating derivatives by central differencing in the y-direction
* `upwind_x(var,v):` a method for calculating derivatives by upwind differencing in the x-direction
* `upwind_y(var,v):` a method for calculating derivatives by upwind differencing in the y-direction
* `timestep():` a method which calcualtes an a timestep length from discretized expressions not yet calculated
* `hydrosolver():` calculates discretized expressions, calculates optimal timestep length, and evolves the system in time
* `add_temperature_perturbation(temperature_peak, x0, y0, spread_x, spread_y):` adds a disturbance in the initial equilibrium condition of temperature. This disturbance is Gaussian, and can be placed with a centre in an arbitrary spot in the xy-space, with an arbitrary intensity `temperature_peak`.
