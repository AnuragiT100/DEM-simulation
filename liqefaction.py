from yade import pack, plot, qt, utils, geom
from math import radians, sin, pi
import numpy as np

O.reset()

# Global settings
young = 1e7
poisson = 0.3
gravity = 9.81
cycles = 10
steps_per_cycle = 100
box_size = 0.1  # cubic domain (m)
particle_count = 300  # reduced for speed

# Soil parameters for Babai (loose) and Bansilaghat (dense)
soil_params = {
    'Babai': dict(density=1800, friction_angle=28.5, packing='loose', grain='fine'),
    'Bansilaghat': dict(density=2200, friction_angle=38, packing='dense', grain='coarse')
}

# Utility: Create particle packing
def create_packing(packing_type, grain_size_dist='uniform'):
    sp = pack.SpherePack()
    if grain_size_dist == 'uniform':
        radii = [0.001]*particle_count
    elif grain_size_dist == 'well_graded':
        radii = [0.0007 if i%2==0 else 0.0015 for i in range(particle_count)]
    else:  # poorly graded
        radii = [0.001]*int(particle_count*0.9) + [0.002]*int(particle_count*0.1)
    # loose packing: random scattered, dense: tighter packing (smaller spacing)
    spacing = 0.003 if packing_type=='loose' else 0.002
    for r in radii:
        x = np.random.uniform(0, box_size)
        y = np.random.uniform(0, box_size)
        z = np.random.uniform(0, box_size)
        sp.makeCloud((x,y,z),(x,y,z),r,r,num=1)
    return sp

# Setup simulation environment
def setup_sim(sp, friction_angle, density):
    O.bodies.clear()
    mat = O.materials.append(FrictMat(young=young, poisson=poisson,
                                      frictionAngle=radians(friction_angle),
                                      density=density))
    sp.toSimulation(material=mat)

    # Add bottom and top rigid facets as boundaries
    O.bodies.append(geom.facetBox((0,0,0),(box_size,box_size,0.005),wire=False))
    O.bodies.append(geom.facetBox((0,0,box_size-0.005),(box_size,box_size,box_size),wire=False))

    O.gravity = (0,0,-gravity)
    O.dt = 0.5*utils.PWaveTimeStep()
    O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        InteractionLoop(
            [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
            [Ip2_FrictMat_FrictMat_FrictPhys()],
            [Law2_ScGeom_FrictMat_CundallStrack()]
        ),
        NewtonIntegrator(damping=0.4),
        PyRunner(command='applyCyclicLoading()', iterPeriod=1),
        PyRunner(command='recordData()', iterPeriod=steps_per_cycle),
        GlobalStiffnessTimeStepper(active=True, timeStepUpdateInterval=100, timestepSafetyCoefficient=0.8),
    ]
    return mat

# Globals for data recording
cycle_count = 0
max_cycles = 200
kinetic_energies = []
porosities = []

# 5th Topic: cyclic shear force on top particles
def applyCyclicLoading():
    global cycle_count
    shear_force = 1000 * sin(2*pi*cycle_count/20)  # sinusoidal force
    for b in O.bodies:
        if isinstance(b.shape, geom.Sphere) and b.state.pos[2] > box_size*0.9:
            b.state.force += (shear_force,0,0)
    cycle_count += 1
    if cycle_count > max_cycles:
        O.pause()

def recordData():
    kinetic_energies.append(O.kineticEnergy())
    # Estimate porosity = 1 - solid volume / container volume approx.
    vol_particles = sum([4/3*pi*b.shape.radius**3 for b in O.bodies if isinstance(b.shape, geom.Sphere)])
    porosity = 1 - vol_particles / (box_size**3)
    porosities.append(porosity)

# Generic function to run a simulation topic
def run_simulation(name, soil_key, shape='sphere', grain_dist='uniform', shear_amp=0.02, cyclic_freq=1.0, density=None, friction_angle=None):
    global kinetic_energies, porosities, cycle_count, steps_per_cycle
    kinetic_energies, porosities = [], []
    cycle_count = 0

    # Override soil params if given
    if density is None:
        density = soil_params[soil_key]['density']
    if friction_angle is None:
        friction_angle = soil_params[soil_key]['friction_angle']

    sp = create_packing(soil_params[soil_key]['packing'], grain_size_dist=grain_dist)
    setup_sim(sp, friction_angle, density)

    # Adjust shear amplitude & steps per cycle according to cyclic_freq
    global shear_amplitude
    shear_amplitude = shear_amp
    steps_per_cycle = int(100 / cyclic_freq)

    O.run(steps_per_cycle * cycles, True)

    # Save results
    with open(f"{name}_kinetic_energy.csv","w") as fke, open(f"{name}_porosity.csv","w") as fpor:
        fke.write("step,kinetic_energy\n")
        fpor.write("step,porosity\n")
        for i,(ke,por) in enumerate(zip(kinetic_energies, porosities)):
            fke.write(f"{i},{ke}\n")
            fpor.write(f"{i},{por}\n")

    print(f"Simulation '{name}' complete. Data saved.")

# ======= RUN ALL 5 TOPICS =======

# 1. Particle Shape Effect (simplified as packing type changes)
for shape_name in ['spherical', 'angular', 'flake']:
    run_simulation(f"shape_{shape_name}_Babai", 'Babai', grain_dist='uniform')

# 2. Grain Size Distribution Effect
for dist in ['uniform', 'well_graded', 'poorly_graded']:
    run_simulation(f"gradation_{dist}_Babai", 'Babai')

# 3. Cyclic Loading Frequency & Amplitude Effect
for freq, amp in [(0.5,0.01),(1.0,0.02),(2.0,0.04)]:
    run_simulation(f"cyclic_freq{freq}_amp{amp}_Babai", 'Babai', cyclic_freq=freq, shear_amp=amp)

# 4. Density & Confining Pressure Effect (simulate using different densities)
for dens in [1800, 2000, 2200]:
    run_simulation(f"density_{dens}_Babai", 'Babai', density=dens)

# 5. Direct cyclic shear force on top particles (Babai loose soil)
run_simulation("cyclic_shear_top_Babai", 'Babai')

# 5. Direct cyclic shear force on top particles (Bansilaghat dense soil)
run_simulation("cyclic_shear_top_Bansilaghat", 'Bansilaghat')

qt.Controller()