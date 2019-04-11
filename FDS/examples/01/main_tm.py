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
app1.initialize(metaData=md)
app2 = apis_and_linking.oofem_api()
app2.initialize(metaData=md, file="input.tm")

app1.meshBCID = 1

filenames = []
nsteps = app1.loadDumpedInfo("FDSASTF.tpf", filenames)
tmpFld = app1.loadASTField(filenames[0])
app2.ASTField.createFromASTField(tmpFld)

time = mupif.Physics.PhysicalQuantities.PhysicalQuantity(0.0, 's')
time_prev = mupif.Physics.PhysicalQuantities.PhysicalQuantity(0.0, 's')
dt = mupif.Physics.PhysicalQuantities.PhysicalQuantity(0.0, 's')

print(nsteps)

app2.regASTField()

for i in range(0, nsteps):
    tmpFld = app1.loadASTField(filenames[i])
    time_prev = time
    time = mupif.Physics.PhysicalQuantities.PhysicalQuantity(tmpFld.time, 's')
    time_step = mupif.TimeStep.TimeStep(time, dt, time + dt, n=i+1)
    app2.setASTField(app1.giveASTFieldValues())
    app2.solveStep(time_step)

app2.terminate()
