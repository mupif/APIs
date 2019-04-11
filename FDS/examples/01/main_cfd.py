import sys
import mupif
sys.path.append('../..')
import apis_and_linking

md = {
    'Execution': {
        'ID': '1',
        'Use_case_ID': '1',
        'Task_ID': '1'
    }
}

app1 = apis_and_linking.fds_api()
app1.initialize(metaData=md, file="input.fds")

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
