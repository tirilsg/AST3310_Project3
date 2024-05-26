# AST3310_Project3

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

The program contains code which uses this class to run simulations and visualize this data by calling the module `FVIs3.py`, also located in the repository. Which simulation is ran is determined by setting the following boolean variables either True or False; 
* `Sanity == True:` checks the sanity of the initialization method by plotting the temperature distibution in the 2D space at different times, t=0 and t=60 seconds.
* `Single == True:` simualtes the behaviours occurring within the space when the system is initialised with a single disturbance in the temperature distribution. The simulation is ran for a duration of 600s, with data saved in a map `"single_disturbances"`. The initial condition of both temperature and density spread is plotted and saved as `"figures/initial_two_single.pdf"`, and the behaviours of temperature, velocity, density and enerfy flux is animated. Lastly, the behaviour of the average of all relevant variables internal energy e, mommentum density rv, energy flux ef, temperature T, mass density $rho$ and speed v is visualized over time. These plots and animations are obtained by setting the perturbation_sim arguments `velocities=True`, `densities=True`, `avg = True`.
* `Small_single == True:` initializes the system with a cooler perturbation in temperature than previously done.  The initial condition of both temperature and density spread is plotted and saved as `"figures/initial_two_single_small.pdf"`.
* `Large_single == True:` creates the same simulation as `Single == True`, with a wider spread of temperature disturbances in the x-direction
* `Five == True:` initializes a system with five disturbances in the initial temperature distributions. The simulation is ran for a duration of 600s, with data saved in a map `"five_disturbances"`. Animations are created for both energy flux and temperature distribuition.
* `Six == True:` initializes a system with five disturbances in the initial temperature distributions. The simulation is ran for a duration of 600s, with data saved in a map `"six_disturbances"`. Animations are created for both energy flux and temperature distribuition.
