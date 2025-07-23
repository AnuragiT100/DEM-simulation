

"""
bearing_capacity_test.py
------------------------
Simple YADE simulation for preliminary plate load test
to estimate bearing capacity trends for borehole soils BH-5 (granular)
and BH-6 (clayey).

Author: Anuragi Thapa
Date-2082-03-16
"""

from yade import pack, plot, qt, utils
from yade import FrictMat, ForceResetter, InsertionSortCollider, InteractionLoop, \
    Ig2_Sphere_Sphere_ScGeom, Ig2_Box_Sphere_ScGeom, Ip2_FrictMat_FrictMat_FrictPhys, \
    Law2_ScGeom_FrictPhys_CundallStrack, NewtonIntegrator, PyRunner, PWaveTimeStep
from math import radians

# --- Material definitions ---
granularSoil = FrictMat(young=1e7, poisson=0.3, frictionAngle=radians(30), density=2600, label='granular')
clayeySoil = FrictMat(young=1e6, poisson=0.4, frictionAngle=radians(20), density=1800, label='clayey')
plateMat = FrictMat(young=2e8, poisson=0.25, frictionAngle=0, density=7850, label='plate')

# --- Choose soil type ---
# For BH-5 (granular), use granularSoil
# For BH-6 (clayey/mixed), use clayeySoil
soilMat = granularSoil  # Change to clayeySoil for BH-6

# Add materials to simulation
O.materials.append(granularSoil)
O.materials.append(clayeySoil)
O.materials.append(plateMat)

# --- Create soil packing ---
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (1, 1, 0.5), rMean=0.02, rRelFuzz=0.3, num=5000)
sp.toSimulation(material=soilMat)

# --- Create rigid plate on top ---
plate = utils.box(center=(0.5, 0.5, 0.55), extents=(0.2, 0.2, 0.01), fixed=True, wire=True, material=plateMat)
plateId = O.bodies.append(plate)

# --- Apply constant downward load (plate pressure) ---
load = -1e4  # Newtons (negative for downward force)
O.forces.setPermF(plateId, (0, 0, load))

# --- Define simulation engines ---
O.engines = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    NewtonIntegrator(gravity=(0, 0, -9.81), damping=0.4),
    PyRunner(iterPeriod=100, command='recordData()'),
]

# --- Set time step ---
O.dt = 0.5 * PWaveTimeStep()

# --- Data recording function ---
def recordData():
    dispZ = O.bodies[plateId].state.displ()[2]  # vertical displacement of plate
    plot.addData(step=O.iter, displacement=-dispZ)  # negative to show settlement as positive

# --- Start simulation ---
O.run(5000, True)  # Run 5000 steps or until user stops

# --- Plot results at the end ---
qt.Controller().show()
plot.saveDataTxt('bearing_capacity_results.txt')
plot.saveFig('bearing_capacity_plot.png')

