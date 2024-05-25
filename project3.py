import numpy as np 
import FVis3 as FVis
import matplotlib.pyplot as plt 
import shutil ; import os

class Convection2D:
    def __init__(self):
        self.rad_sun = 6.96*10**8         # units meters, solar radius
        self.mass_sun = 1.989*10**30      # units of kg, solar mass
        
        #definitions of constants and units
        self.u_val = 1.660539*10**(-27)   # unit kg, atomic mass unit u convertion
        self.kb = 1.3806*10**(-23)        # unit m^2kg/(s^2K), Boltzmann's constant 
        self.G = 6.6742*10**(-11)         # unit of Nm^2/kg^2, gravitational constant for the sun
        self.c = 2.9979*10**8             # unit of m/s, speed of light
        
        self.g = self.G * self.mass_sun / (self.rad_sun**2)
        self.mean_molecular = 0.61        # mean molecular weight
        self.gradient = 1/2               # gradient redefined 
        self.pressure = 1.8*10**4         # unit Pa, solar photosphere pressure 
        self.temperature = 5775           # unit K, solar photosphere temperature
        self.gamma = 5/3                  # gamma for an ideal gas
        
        #box-defining variables
        self.X=12e6            ; self.Y=4e6
        self.nx=300            ; self.ny=100
        self.dx=self.X/self.nx ; self.dy=self.Y/self.ny
        self.p=0.1 
        
        #arrays to be filled with data
        self.T = np.empty((self.ny,self.nx), dtype=float)
        self.P = np.empty((self.ny,self.nx), dtype=float)
        self.u = np.zeros((self.ny,self.nx), dtype=float)
        self.w = np.zeros((self.ny,self.nx), dtype=float)
        self.rho = np.empty((self.ny,self.nx), dtype=float)
        self.e = np.empty((self.ny,self.nx), dtype=float)
        
        #arrays to contain discretized expressions
        self.drho_dt= np.empty((self.ny, self.nx), dtype=float)
        self.drhou_dt= np.empty((self.ny, self.nx), dtype=float)
        self.drhow_dt = np.empty((self.ny, self.nx), dtype=float)
        self.de_dt= np.empty((self.ny, self.nx), dtype=float)
        
        #arbitrary gaussian temperature disturbance
        self.arb_pert = np.zeros((self.ny, self.nx))
    
    #method that initializes the system 
    def initialise(self):
        self.y   = np.linspace(0, self.Y, self.ny)
        for i in range(self.ny):
            self.T[i , :] = (self.g * self.gradient * self.u_val * self.mean_molecular / self.kb * (self.Y - self.y[i])) +self.temperature 
        self.P   = self.pressure*(self.T / self.temperature)**(1 / self.gradient)
        self.T  += self.arb_pert
        self.e   = self.P / (self.gamma - 1)
        self.rho = self.e * (self.gamma - 1) * self.mean_molecular * self.u_val / (self.kb * self.T)
        
    #implementing boundaries for the behaviours of the gas at the edge of the 2D-space  
    def boundary_conditions(self):
        self.w[0 , :]    = 0
        self.w[-1 , :]   = 0
        self.u[0 , :]    = (4 * self.u[1 , :] - self.u[2 , :]) / 3
        self.u[-1 , :]   = (4 * self.u[-2 , :] - self.u[-3 , :]) / 3
        self.e[0 , :]    = (4 * self.e[1 , :] - self.e[2 , :]) / (3 - 2*self.dy * ((self.g*self.mean_molecular*self.u_val) / (self.kb*self.T[0 , :])))
        self.e[-1 , :]   = (4 * self.e[-2 , :] - self.e[-3 , :]) / (3 + 2*self.dy * ((self.g*self.mean_molecular*self.u_val) / (self.kb*self.T[-1 , :])))
        self.rho[0 , :]  = (self.gamma - 1) * self.mean_molecular  * self.u_val * self.e[0 , :] / (self.kb * self.T[0 , :])
        self.rho[-1 , :] = (self.gamma - 1) * self.mean_molecular * self.u_val * self.e[-1 , :] / (self.kb * self.T[-1 , :])
        
    def central_x(self, var):
        previous = np.roll(var, -1, 1)
        next= np.roll(var, 1, 1)
        return ((previous - next) / (2 * self.dx))

    def central_y(self, var):
        previous = np.roll(var, -1, 0)
        next= np.roll(var, 1, 0)
        return ((previous - next) / (2 * self.dy))

    def upwind_x(self, var, v):
        next = np.roll(var, 1, 1)
        previous = np.roll(var, -1, 1)
        return (np.where(v < 0, (previous - var) / self.dx, (var-next) / self.dx))

    def upwind_y(self, var, v):
        next = np.roll(var, 1, 0)
        previous= np.roll(var, -1, 0)
        return (np.where(v < 0, (previous-var)/self.dy, (var-next) / self.dy))

    # defining time step size
    def timestep(self):
        p = 0.1
        r_u       = np.nanmax(np.abs(self.u/self.dx))
        r_w       = np.nanmax(np.abs(self.w/self.dy))
        r_e       = np.nanmax(np.abs(self.de_dt/self.e))
        r_rho     = np.nanmax(np.abs(self.drho_dt/self.rho))
        delta = np.nanmax([r_rho, r_u, r_w, r_e])
        if delta == 0: #to avoid division by zero, the 0 is replaced by 1
            delta = 1
        self.dt = p / delta

    # function that defines discetized behaviours and evolves system in time
    def hydro_solver(self):
        ######################### Continuity #########################
        dp_dx          = self.central_x(self.P)
        dp_dy          = self.central_y(self.P)
        du_dx          = self.central_x(self.u)
        dw_dy          = self.central_y(self.w)
        drho_dx        = self.upwind_x(self.rho, self.u)
        drho_dy        = self.upwind_y(self.rho, self.w)
        self.drho_dt   = -self.rho* (du_dx + dw_dy) - self.u * drho_dx - self.w * drho_dy

        #################### Horizontal  Momentum ####################
        self.rhou      = self.rho * self.u
        du_dx          = self.upwind_x(self.u, self.u)
        drhou_dx       = self.upwind_x(self.rhou,self.u)
        drhou_dy       = self.upwind_y(self.rhou,self.w)
        self.drhou_dt  = -(self.rhou) * (du_dx + dw_dy) - self.u * drhou_dx - self.w * drhou_dy - dp_dx

        ##################### Vertical  Momentum #####################
        self.rhow      = self.rho * self.w
        du_dx          = self.central_x(self.u)
        dw_dy          = self.upwind_y(self.w, self.w)
        drhow_dx       = self.upwind_x(self.rhow,self.u)
        drhow_dy       = self.upwind_y(self.rhow,self.w)
        self.drhow_dt = -(self.rhow) * (du_dx + dw_dy) - self.u * drhow_dx - self.w * drhow_dy - dp_dy - self.rho * self.g

        ########################### Energy ###########################
        du_dx          = self.central_x(self.u)
        dw_dy          = self.central_y(self.w)
        de_dx          = self.upwind_x(self.e, self.u)
        de_dy          = self.upwind_y(self.e, self.w)
        self.de_dt     = -self.u * de_dx - self.w * de_dy - (self.P+ self.e) * (du_dx + dw_dy)
    
        ##################### Evolution in Time #####################
        self.timestep() #defining the size of the time step 
        self.e[:] = self.e + self.de_dt * self.dt
        self.rho[:] = self.rho + self.drho_dt * self.dt
        self.u[:] = (self.rho * self.u + self.drhou_dt * self.dt) / self.rho
        self.w[:] = (self.rho * self.w + self.drhow_dt * self.dt) / self.rho
        #redefine boundaries before calculating pressure and temperature
        self.boundary_conditions() 
        self.P[:] = (self.gamma - 1) * self.e 
        self.T[:] = self.P * self.u_val *self.mean_molecular / (self.rho * self.kb)
        return self.dt
    
    # temperature disturbance definition
    def add_temperature_perturbation(self, temperature_peak, x0, y0, spread_x, spread_y):
        x_array = np.linspace(0, self.X, self.nx)
        y_array = np.linspace(0, self.Y, self.ny)
        x, y = np.meshgrid(x_array, y_array)
        #the disturbance in the temperature is saved as a sum, so that multiple perturbations can be added to a system if needed
        self.arb_pert += temperature_peak * np.exp(-((x - x0)**2/(2*spread_x**2) + (y-y0)**2/(2*spread_y**2)))  

def plot_initial_conditions(instance, name, variable_name):
    if variable_name=="two":
        fig, axs = plt.subplots(1, 2, figsize=(15, 5))
        axs[0].set_xlabel("x")
        axs[0].set_ylabel("y")
        axs[0].set_title(f"Temperature")
        fig.colorbar(axs[0].pcolormesh(instance.T[:], cmap="plasma", shading="auto"), ax=axs[0])
        axs[1].set_xlabel("x")
        axs[1].set_ylabel("y")
        axs[1].set_title(f"Density")
        fig.colorbar(axs[1].pcolormesh(instance.rho[:], cmap="plasma", shading="auto"), ax=axs[1])
        fig.tight_layout()
    else:
        fig, axs = plt.subplots(figsize=(8, 6))
        axs.set_xlabel("x")
        axs.set_ylabel("y")
        axs.set_title(f"{variable_name}")
        if variable_name=="Temperature":
            fig.colorbar(axs.pcolormesh(instance.T[:], cmap="plasma", shading="auto"))
        if variable_name=="Pressure":
            fig.colorbar(axs.pcolormesh(instance.P[:], cmap="plasma", shading="auto"))
        if variable_name=="Density":
            fig.colorbar(axs.pcolormesh(instance.rho[:], cmap="plasma", shading="auto"))
        if variable_name=="Internal Energy":
            fig.colorbar(axs.pcolormesh(instance.e[:], cmap="plasma", shading="auto"))
        fig.tight_layout()
    plt.savefig(f"{name}")
    plt.show()

def perturbation_sim(instance, shot_n, time, foldername, quiverscale, velocities=False, densities=False, avg = False):
    fluid = FVis.FluidVisualiser()
    units = {"Lx":"Mm","Lz":"Mm"}
    snapshots =  list(np.arange(0,time+1,shot_n))
    fluid.save_data(time, instance.hydro_solver, rho=instance.rho, T=instance.T, u=instance.u,  w=instance.w,P=instance.P,e=instance.e,sim_fps=30,folder=foldername)
    fluid.animate_2D('T', height=4.6, quiverscale=quiverscale, snapshots = snapshots, units=units, cmap="plasma", showQuiver=True, folder=f"{foldername}", anim_fps=0.15)
    fluid.animate_2D('T', height=4.6, quiverscale=0.35, units=units, cmap="plasma", showQuiver=True, folder=f"{foldername}", anim_fps=0.15)
    shutil.copytree(foldername, f"{foldername}ef") #copying the map and renaming in order to take snapshots with different name indicating energy flux
    fluid.animate_energyflux(folder=f"{foldername}ef" ,height=4.6, snapshots = snapshots, units=units, anim_fps=0.15)
    fluid.animate_energyflux(folder=f"{foldername}ef" ,height=4.6, units=units, anim_fps=0.15)
    if densities == True:
        shutil.copytree(foldername, f"{foldername}rho") #copying the map and renaming in order to take snapshots with different name indicating density
        fluid.animate_2D('rho', height=4.6, quiverscale=quiverscale, snapshots = snapshots, units=units, cmap="plasma", showQuiver=True, folder=f"{foldername}rho", anim_fps=0.15)
        fluid.animate_2D('rho', height=4.6, quiverscale=quiverscale, units=units, cmap="plasma", showQuiver=True, folder=f"{foldername}rho", anim_fps=0.15)
    if velocities == True:
        shutil.copytree(foldername, f"{foldername}v") #copying the map and renaming in order to take snapshots with different name indicating velocity
        fluid.animate_2D('v', height=4.6, quiverscale=quiverscale, snapshots = snapshots, units=units, cmap="plasma",showQuiver=True, folder=f"{foldername}v", anim_fps=0.15)
        fluid.animate_2D('v', height=4.6, quiverscale=quiverscale, units=units, cmap="plasma",showQuiver=True, folder=f"{foldername}v", anim_fps=0.15)
    if avg == True:
        #the method plot_avg does not allow for automatic saving, so each figure must be saved manually
        fluid.plot_avg('rv', relative=True, showTrendline=True, units=units, folder=f"{foldername}")
        fluid.plot_avg('v', relative=True, showTrendline=True, units=units, folder=f"{foldername}")
        fluid.plot_avg('ev', relative=True, showTrendline=True, units=units, folder=f"{foldername}")
        fluid.plot_avg('T', relative=True, showTrendline=True, units=units, folder=f"{foldername}")
        fluid.plot_avg('rho', relative=True, showTrendline=True, units=units, folder=f"{foldername}")
        fluid.plot_avg('e', relative=True, showTrendline=True, units=units, folder=f"{foldername}")


#sanity check for the simulation over time
Sanity = True
if Sanity ==True:
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    instance = Convection2D()
    instance.initialise()
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    axs[0].set_title(f"Temperature, t=0s")
    fig.colorbar(axs[0].pcolormesh(instance.T[:], cmap="plasma", shading="auto"), ax=axs[0])
    fluid = FVis.FluidVisualiser()
    fluid.save_data(60, instance.hydro_solver, rho=instance.rho, T=instance.T, u=instance.u,  w=instance.w,P=instance.P,e=instance.e,sim_fps=30,folder="no_disturbances")
    axs[1].set_xlabel("x")
    axs[1].set_ylabel("y")
    axs[1].set_title(f"Temperature, t=60s")
    fig.colorbar(axs[1].pcolormesh(instance.T[:], cmap="plasma", shading="auto"), ax=axs[1])
    fig.tight_layout()
    plt.savefig(f"figures/sanity_check.pdf")
    plt.show()
 
#simulation of a single temperature disturbance over 600s
Single =True
if Single == True:
    instance = Convection2D()
    instance.add_temperature_perturbation(temperature_peak=60000, x0=6e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.initialise()
    plot_initial_conditions(instance, "figures/initial_two_single.pdf","two")
    perturbation_sim(instance, 30, 600, foldername = "single_disturbances", quiverscale=0.35, velocities=True, densities=True, avg=True)

#simulation of five disturbances over 600s
Five=True
if Five == True:
    instance = Convection2D()
    instance.add_temperature_perturbation(temperature_peak=60000, x0=2e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.add_temperature_perturbation(temperature_peak=40000, x0=4e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.add_temperature_perturbation(temperature_peak=60000, x0=6e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.add_temperature_perturbation(temperature_peak=40000, x0=8e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.add_temperature_perturbation(temperature_peak=60000, x0=10e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.initialise()
    plot_initial_conditions(instance, "figures/initial_two_five.pdf","two")
    perturbation_sim(instance, 30, 600, foldername = "five_disturbances", quiverscale=0.35)
    
#initialisation of a single, cooler perturbation 
Small_single = True
if Small_single == True:
    instance = Convection2D()
    instance.add_temperature_perturbation(temperature_peak=8000, x0=6e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.initialise()
    plot_initial_conditions(instance, "figures/initial_two_single_small.pdf","two")
    
#simulation of a single, wider, disturbance over 200s
Large_single = False
if Large_single == True:
    instance = Convection2D()
    instance.add_temperature_perturbation(temperature_peak=60000, x0=6e6, y0=0, spread_x=1.5e6, spread_y=4e6)
    instance.initialise()
    plot_initial_conditions(instance, "figures/initial_two_single_large.pdf","two")
    perturbation_sim(instance, 10, 200, foldername = "single_disturbances_large",quiverscale=0.6)

#simulation of six disturbances over 600s
Six = False
if Six == True:
    instance = Convection2D()
    instance.add_temperature_perturbation(temperature_peak=60000, x0=2e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.add_temperature_perturbation(temperature_peak=40000, x0=4e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.add_temperature_perturbation(temperature_peak=60000, x0=5.5e6, y0=0, spread_x=5e5, spread_y=5e6)
    instance.add_temperature_perturbation(temperature_peak=60000, x0=6.5e6, y0=0, spread_x=5e5, spread_y=5e6)
    instance.add_temperature_perturbation(temperature_peak=40000, x0=8e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.add_temperature_perturbation(temperature_peak=60000, x0=10e6, y0=0, spread_x=5e5, spread_y=3e6)
    instance.initialise()
    plot_initial_conditions(instance, "figures/initial_two_six.pdf","two")
    perturbation_sim(instance, 30, 600, foldername = "six_disturbances",quiverscale=0.35)
