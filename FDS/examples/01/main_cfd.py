import sys
import mupif
sys.path.append('../..')
import apis_and_linking as fol

app1 = fol.fds_api("input.fds")
app1.initialize()

app1.meshBCID = 1

app1.loadMeshes()
app1.loadASTMesh()

time = mupif.Physics.PhysicalQuantities.PhysicalQuantity(0.0, 's')
targetTime = mupif.Physics.PhysicalQuantities.PhysicalQuantity(20.5, 's')
dt = mupif.Physics.PhysicalQuantities.PhysicalQuantity(2.0, 's')

step = 0
while time < targetTime:
    step += 1
    time = time + dt
    time_step = mupif.TimeStep.TimeStep(time, dt, time + dt, n=step)
    app1.solveStep(time_step)
    app1.exportMeshes()
    app1.loadASTemps()
    app1.dumpASTField()

app1.terminate()
