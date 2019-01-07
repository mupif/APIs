from __future__ import print_function, division
from mupif import *
import numpy as np
from ctypes import *
# from numpy.ctypeslib import ndpointer
import os

oofemImported = False
liboofem = None

c_double_p = POINTER(c_double)
c_int_p = POINTER(c_int)
errtol = 0.00001

time = 0.0
targetTime = 0.0
fdsTStep = 0
oofemTStep = 0


class BrickShape:
    def __init__(self, dimensions):
        self.corners = []
        for i in range(0, 6):
            self.corners.append(dimensions[i])


class ClassCentroid:
    def __init__(self, centroids):
        self.centroids = centroids

    def isInsideBox(self, mesh_v):
        b1 = self.centroids[0] > mesh_v.corners[0]
        b2 = self.centroids[0] < mesh_v.corners[3]
        b3 = self.centroids[1] > mesh_v.corners[1]
        b4 = self.centroids[1] < mesh_v.corners[4]
        b5 = self.centroids[2] > mesh_v.corners[2]
        b6 = self.centroids[2] < mesh_v.corners[5]
        if b1 and b2 and b3 and b4 and b5 and b6:
            return 1
        return 0


class TempField:
    def __init__(self, numNodes=0, numCells=0):
        global oofemImported
        global liboofem
        if not oofemImported:
            oofemImported = True
            import liboofem as liboofem

        self.uniform = None  # 1=uniform, 0=unstructured, -1=field of points without cells
        self.TField = None
        self.nnodes = None
        self.ncells = None
        self.meshSize = [0, 0, 0]
        self.domainSize = [0., 0., 0., 0., 0., 0.]

    def setType(self, uniform, numNodes=0, numCells=0):
        if uniform > 0:
            self.uniform = 1
        else:
            if uniform == 0:
                self.uniform = 0
            else:
                self.uniform = -1

        if self.uniform == 1:
            self.TField = liboofem.UniformGridField()
        else:
            if uniform == 0:
                self.TField = liboofem.UnstructuredGridField(numNodes, numCells)
                self.nnodes = numNodes
                self.ncells = numCells
            else:
                self.TField = liboofem.UnstructuredGridField(numNodes, 0)
                # or change the type of field when implemented
                self.nnodes = numNodes
                self.ncells = 0

    def setValues(self, valField):
        if self.uniform == 1:
            oField = liboofem.FloatArray(len(valField))
            for i in range(0, len(valField)):
                oField[i] = float(valField[i])
            self.TField.setValues(oField)  # it is a different "setValues"...is is an oofem function
        else:
            if self.uniform == 0 or self.uniform == -1:
                for i in range(0, self.nnodes):
                    varVal = liboofem.FloatArray(1)
                    varVal[0] = float(valField[i])
                    self.TField.setVertexValue(i, varVal)

    def createFromASTField(self, someField):
        TotMeshNNodes = len(someField.mesh.vertexList)
        self.setType(-1, TotMeshNNodes)
        for i in range(0, TotMeshNNodes):
            coo = liboofem.FloatArray(3)
            coo[0] = someField.mesh.vertexList[i].coords[0]
            coo[1] = someField.mesh.vertexList[i].coords[1]
            coo[2] = someField.mesh.vertexList[i].coords[2]
            self.TField.addVertex(i, coo)

    def createFromTempField(self, mshID, someField, cellType=1):
        if mshID == 0:
            TotMeshNNodes = len(someField.mesh.vertexList)
            TotMeshNCells = len(someField.mesh.cellList)
            if cellType == 1:
                self.setType(0, TotMeshNNodes, TotMeshNCells)
            if cellType == 0:
                self.setType(0, TotMeshNNodes, TotMeshNCells * 6)
            for i in range(0, TotMeshNNodes):
                coo = liboofem.FloatArray(3)
                coo[0] = someField.mesh.vertexList[i].coords[0]
                coo[1] = someField.mesh.vertexList[i].coords[1]
                coo[2] = someField.mesh.vertexList[i].coords[2]
                self.TField.addVertex(i, coo)
            for i in range(0, TotMeshNCells):
                cvlp = someField.mesh.cellList[i].vertices
                if cellType == 1:
                    verts = liboofem.IntArray(8)
                    for j in range(0, 8):
                        verts[j] = cvlp[j]
                    self.TField.addCell(i, liboofem.Element_Geometry_Type.EGT_hexa_1, verts)

                if cellType == 0:
                    verts = liboofem.IntArray(4)
                    verts[0] = cvlp[0]
                    verts[1] = cvlp[3]
                    verts[2] = cvlp[4]
                    verts[3] = cvlp[1]
                    self.TField.addCell(6 * i, liboofem.Element_Geometry_Type.EGT_tetra_1, verts)
                    verts = liboofem.IntArray(4)
                    verts[0] = cvlp[4]
                    verts[1] = cvlp[1]
                    verts[2] = cvlp[3]
                    verts[3] = cvlp[2]
                    self.TField.addCell(6 * i + 1, liboofem.Element_Geometry_Type.EGT_tetra_1, verts)
                    verts = liboofem.IntArray(4)
                    verts[0] = cvlp[4]
                    verts[1] = cvlp[1]
                    verts[2] = cvlp[2]
                    verts[3] = cvlp[5]
                    self.TField.addCell(6 * i + 2, liboofem.Element_Geometry_Type.EGT_tetra_1, verts)
                    verts = liboofem.IntArray(4)
                    verts[0] = cvlp[4]
                    verts[1] = cvlp[3]
                    verts[2] = cvlp[7]
                    verts[3] = cvlp[6]
                    self.TField.addCell(6 * i + 3, liboofem.Element_Geometry_Type.EGT_tetra_1, verts)
                    verts = liboofem.IntArray(4)
                    verts[0] = cvlp[4]
                    verts[1] = cvlp[3]
                    verts[2] = cvlp[6]
                    verts[3] = cvlp[2]
                    self.TField.addCell(6 * i + 4, liboofem.Element_Geometry_Type.EGT_tetra_1, verts)
                    verts = liboofem.IntArray(4)
                    verts[0] = cvlp[4]
                    verts[1] = cvlp[2]
                    verts[2] = cvlp[6]
                    verts[3] = cvlp[5]
                    self.TField.addCell(6 * i + 5, liboofem.Element_Geometry_Type.EGT_tetra_1, verts)
        else:
            self.setType(1)
            # Zjistim minmax a deleni numericky. Prasarna, ale univerzalni.
            locerrtol = 0.01
            self.domainSize = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            self.meshSize = [0, 0, 0]
            a = someField.mesh.vertexList[0].coords
            b = someField.mesh.vertexList[len(someField.mesh.vertexList) - 1].coords
            self.domainSize = [a[0], a[1], a[2], b[0], b[1], b[2]]

            for i in range(0, len(someField.mesh.vertexList)):
                for j in range(0, 3):
                    inLine = 1
                    for k in range(0, 3):
                        if j != k:
                            if (someField.mesh.vertexList[i].coords[k] > a[k] + locerrtol or
                                    someField.mesh.vertexList[i].coords[k] < a[k] - locerrtol):
                                inLine = 0
                    if inLine:
                        self.meshSize[j] = self.meshSize[j] + 1

            boxSize = self.domainSize
            meshDiv = self.meshSize
            self.TField.setGeometry(lo=(boxSize[0], boxSize[1], boxSize[2]), hi=(boxSize[3], boxSize[4], boxSize[5]),
                                    div=(meshDiv[0] - 1, meshDiv[1] - 1, meshDiv[2] - 1))
            self.nnodes = meshDiv[0] * meshDiv[1] * meshDiv[2]
            self.nnodes = (meshDiv[0] - 1) * (meshDiv[1] - 1) * (meshDiv[2] - 1)


class fds_runall:
    def __init__(self):
        self.libfds = CDLL("/home/stanislav/shared/libfds_v2.so")
        self.libfds.run_all()


# ##########################################################################################################
# ####     ##     #####     ######   ##   ##   #############################################################
# ####  #####  ##  ###   #########   ##   ##   #############################################################
# ####     ##  ###  ###    #######  ###  ###  ##############################################################
# ####  #####  ##  ######   ################################################################################
# ####  #####     ####     ######  ###  ###  ###############################################################
# ##########################################################################################################
class fds_api(Application.Application):
    def __init__(self, filename=""):
        Application.Application.__init__(self)
        self.nMeshes = 0
        self.libfds = None
        if filename != "":
            self.libfds = CDLL("/home/stanislav/shared/libfds_v2.so")  #
            print("FDS reading input file '%s'" % (filename))
            self.libfds.initialize.argtypes = [c_char_p]
            # c_filename=c_char_p(b filename)
            c_filename = bytes(filename, encoding='utf-8')
            self.libfds.initialize(c_filename)
            self.nMeshes = self.libfds.give_nmeshes()
            print('Mesh count: ', self.nMeshes)
        # ID of the mesh, used for the boundary condition.
        self.meshBCID = 1
        self.FDSTempFields = []
        self.FDSASTField = None
        self.nVTU = 0
        self.nVTU_last = 0
        self.VTUoutPartM = False
        self.VTUoutAST = False
        self.times = []
        self.fds_dtime = c_double(0)
        self.fds_time = c_double(0)
        self.step = 0
        self.errorCheck = 0
        self.errorCheckErrTol = 0.00001
        self.errorCheckMeshes = []
        self.errorCheckSteps = []
        self.errorCheckNodes = []
        self.errorCheckValues = []
        self.Meshes = []
        self.ASTMesh = None

    def errorCheckLoad(self, filename="FDSErrorCheck.txt"):
        if os.path.exists(filename):
            print("loading FDSErrorCheck...")
            text_file = open(filename, "r")
            textLine = text_file.readline()
            while len(textLine) > 1:
                splitted = textLine.split()
                if len(splitted) == 4:
                    self.errorCheckMeshes.append(int(splitted[0]))
                    self.errorCheckSteps.append(int(splitted[1]))
                    self.errorCheckNodes.append(int(splitted[2]))
                    self.errorCheckValues.append(float(splitted[3]))
                textLine = text_file.readline()
            text_file.close()
            if len(self.errorCheckValues) > 0:
                self.errorCheck = 1
            text_file = open("errorCheckLog.txt", "a")
            text_file.close()
        else:
            print("FDSErrorCheck file not found")

    def errorCheckRun(self, stepNumber):
        if self.errorCheck:
            for i in range(0, len(self.errorCheckValues)):
                if stepNumber == self.errorCheckSteps[i]:
                    foundValue = float(self.FDSTempFields[self.errorCheckMeshes[i] - 1].value[self.errorCheckNodes[i]])
                    print(" CheckVal %lf .EQ. %lf" % (foundValue, self.errorCheckValues[i]))
                    if self.errorCheckValues[i] == 0.0:
                        print("Zero cannot be checked.")
                    else:
                        computedError = (foundValue - self.errorCheckValues[i]) / self.errorCheckValues[i]
                        if abs(computedError) > self.errorCheckErrTol:
                            if computedError > 0:
                                errorSign = "+"
                            else:
                                errorSign = "-"
                            print("FAIL Error = %s %.2lf %%" % (errorSign, abs(computedError * 100)))
                            # print error to file errorCheckLog.txt
                            text_file = open("errorCheckLog.txt", "a")
                            text_file.write(
                                "-------------------------\nFDS Error = %s %.2lf %%\n   Should be %lf but is %lf\n   Correct form for update: %d %d %d %lf\n\n" % (
                                errorSign, abs(computedError * 100), self.errorCheckValues[i], foundValue,
                                self.errorCheckMeshes[i], self.errorCheckSteps[i], self.errorCheckNodes[i], foundValue))
                            text_file.close()
                        else:
                            print("FDSErrorCheckStepIsOK")

    def getMeshGridProperties(self, meshID, DomainSize, Meshsize):
        ax = c_double(0)
        ay = c_double(0)
        az = c_double(0)
        bx = c_double(0)
        by = c_double(0)
        bz = c_double(0)
        mx = c_int(0)
        my = c_int(0)
        mz = c_int(0)
        paramMeshId = c_int(meshID)
        self.libfds.meshSize(byref(paramMeshId), byref(ax), byref(ay), byref(az), byref(bx), byref(by), byref(bz),
                             byref(mx), byref(my), byref(mz))
        DomainSize[0] = ax.value
        DomainSize[1] = ay.value
        DomainSize[2] = az.value
        DomainSize[3] = bx.value
        DomainSize[4] = by.value
        DomainSize[5] = bz.value
        Meshsize[0] = mx.value
        Meshsize[1] = my.value
        Meshsize[2] = mz.value

    def loadMeshes(self, boundaryTempField=None):
        nnodes = [0, 0, 0]
        dimensions = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        oneCellSizes = []
        meshCellCounts = []
        meshCellCountsBefore = []
        meshCorners = []
        self.meshNNodes = [[0 for x in range(3)] for y in range(self.nMeshes)]  # checkit
        for i in range(0, self.nMeshes):
            newmesh = Mesh.UnstructuredMesh()
            self.getMeshGridProperties(i + 1, dimensions, nnodes)
            meshCorners.append(BrickShape(dimensions))
            self.meshNNodes[i][0] = nnodes[0]
            self.meshNNodes[i][1] = nnodes[1]
            self.meshNNodes[i][2] = nnodes[2]
            nodeid = 0
            ind = [0, 0, 0]
            coords = [0.0, 0.0, 0.0]
            for ind[0] in range(0, nnodes[0]):
                for ind[1] in range(0, nnodes[1]):
                    for ind[2] in range(0, nnodes[2]):
                        for c in range(0, 3):
                            coords[c] = dimensions[c] + ind[c] * (dimensions[c + 3] - dimensions[c]) / (nnodes[c] - 1)
                        newmesh.vertexList.append(Vertex.Vertex(nodeid, nodeid, [coords[0], coords[1], coords[2]]))
                        nodeid = nodeid + 1
            cellid = 0
            ind = [0, 0, 0]
            n = [0, 0, 0, 0, 0, 0, 0, 0]
            yzNodes = nnodes[1] * nnodes[2]
            for ind[0] in range(0, nnodes[0] - 1):
                for ind[1] in range(0, nnodes[1] - 1):
                    for ind[2] in range(0, nnodes[2] - 1):
                        n[0] = 0 + yzNodes * ind[0] + nnodes[2] * ind[1] + ind[2]
                        n[1] = 1 + yzNodes * ind[0] + nnodes[2] * ind[1] + ind[2]
                        n[2] = 1 + yzNodes * ind[0] + nnodes[2] * ind[1] + ind[2] + yzNodes
                        n[3] = 0 + yzNodes * ind[0] + nnodes[2] * ind[1] + ind[2] + yzNodes
                        n[4] = 0 + yzNodes * ind[0] + nnodes[2] * ind[1] + ind[2] + nnodes[2]
                        n[5] = 1 + yzNodes * ind[0] + nnodes[2] * ind[1] + ind[2] + nnodes[2]
                        n[6] = 1 + yzNodes * ind[0] + nnodes[2] * ind[1] + ind[2] + nnodes[2] + yzNodes
                        n[7] = 0 + yzNodes * ind[0] + nnodes[2] * ind[1] + ind[2] + nnodes[2] + yzNodes
                        newmesh.cellList.append(Cell.Brick_3d_lin(newmesh, cellid, cellid, n))
                        cellid = cellid + 1
            meshCellCounts.append(cellid)
            if i > 0:
                meshCellCountsBefore.append(meshCellCountsBefore[i - 1] + meshCellCounts[i - 1])
            else:
                meshCellCountsBefore.append(0)
            velikost = 1.0
            for c in range(0, 3):
                velikost = velikost * (dimensions[c + 3] - dimensions[c]) / (nnodes[c] - 1)
            oneCellSizes.append(velikost)
            self.Meshes.append(newmesh)
            newTempfield = Field.Field(newmesh, FieldID.FID_Temperature, ValueType.ValueType.Scalar, 'C', 0.0)
            self.FDSTempFields.append(newTempfield)

        if boundaryTempField != None:
            boundaryTempField.createFromTempField(self.meshBCID, self.FDSTempFields[self.meshBCID - 1])

    def loadTempsOnMeshes(self):
        tempVal = c_double(0.0)
        for i in range(0, self.nMeshes):
            paramMeshId = c_int(i + 1)
            j = 0
            for mi in range(0, self.meshNNodes[i][0]):
                for mj in range(0, self.meshNNodes[i][1]):
                    for mk in range(0, self.meshNNodes[i][2]):
                        pmi = c_int(mi + 1)
                        pmj = c_int(mj + 1)
                        pmk = c_int(mk + 1)
                        self.libfds.give_teplotu(byref(tempVal), byref(pmi), byref(pmj), byref(pmk), byref(paramMeshId))
                        self.FDSTempFields[i].value[j] = tempVal.value
                        j = j + 1
            self.FDSTempFields[i].time = self.fds_time.value

        if self.VTUoutPartM:
            if self.nVTU_last == self.nVTU:
                self.nVTU = self.nVTU + 1
            for im in range(0, self.nMeshes):
                self.FDSTempFields[im].toVTK3('particularTempField%d_%d.vtu' % (im + 1, self.nVTU))
                pvdFile = open("particularTempField%d.pvd" % (im + 1), "w")
                pvdFile.write("<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n")
                for i in range(0, self.nVTU):
                    pvdFile.write(
                        "<DataSet timestep=\"%lf\" group=\"\" part=\"\" file=\"particularTempField%d_%d.vtu\"/>\n" % (
                        self.times[i], im + 1, i + 1))
                pvdFile.write("</Collection>\n</VTKFile>")
                pvdFile.close()

    def exportMeshes(self):
        self.times.append(self.fds_time.value)
        self.loadTempsOnMeshes()

    def loadASTMesh(self, boundaryTempField=None):
        filename = "./ASTPointsCoords.txt"
        if os.path.exists(filename):
            nodeid = 0
            print("loading AST mesh points...")
            newmesh = Mesh.UnstructuredMesh()
            text_file = open(filename, "r")
            textLine = text_file.readline()
            while len(textLine) > 1:
                splitted = textLine.split()
                if len(splitted) == 3:
                    nodeid = nodeid + 1
                    newmesh.vertexList.append(
                        Vertex.Vertex(nodeid, nodeid, [float(splitted[0]), float(splitted[1]), float(splitted[2])]))
                textLine = text_file.readline()
            text_file.close()
            self.ASTMesh = newmesh
            print("%d ASTPoints" % nodeid)
            self.FDSASTField = Field.Field(self.ASTMesh, FieldID.FID_Temperature, ValueType.ValueType.Scalar, 'C', 0.0)

            if boundaryTempField != None:
                boundaryTempField.createFromASTField(self.FDSASTField)
        else:
            print("ASTPointsCoords.txt file not found")

    def loadASTemps(self):
        for i in range(0, len(self.ASTMesh.vertexList)):
            self.FDSASTField.value[i] = self.giveASTG(i + 1)
        self.FDSASTField.time = self.fds_time.value

        if self.VTUoutAST:
            if self.nVTU_last == self.nVTU:
                self.nVTU = self.nVTU + 1

            self.FDSASTField.toVTK3('AST_field_%d.vtu' % self.nVTU)
            pvdFile = open("AST_field.pvd", "w")
            pvdFile.write("<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n")
            for i in range(0, self.nVTU):
                pvdFile.write(
                    "<DataSet timestep=\"%lf\" group=\"\" part=\"\" file=\"AST_field_%d.vtu\"/>\n" % (self.times[i], i + 1))
            pvdFile.write("</Collection>\n</VTKFile>")
            pvdFile.close()

    def solveStep(self):
        self.nVTU_last = self.nVTU
        self.step = self.step + 1
        self.libfds.solve_step()
        self.libfds.get_t(byref(self.fds_time))
        self.libfds.get_dt(byref(self.fds_dtime))

    def getTime(self):
        return self.fds_time.value

    def complete(self):
        self.libfds.complete()

    def dumpTempFields(self, filename="FDSTempF"):
        if self.nVTU_last == self.nVTU:
            self.nVTU = self.nVTU + 1
        for im in range(0, self.nMeshes):
            if im == self.meshBCID - 1:
                self.FDSTempFields[im].dumpToLocalFile("%s%d_s%d.pf" % (filename, im, self.nVTU), 0)
                text_file = open("%s%d.tpf" % (filename, im), "w")
                text_file.write("%d\n" % self.nVTU)
                for i in range(0, self.nVTU):
                    text_file.write("%s%d_s%d.pf\n" % (filename, im, i + 1))
                text_file.close()

    def loadTempField(self, filename):
        loadedField = Field.Field.loadFromLocalFile(filename)
        if self.meshBCID > len(self.FDSTempFields) - 1:
            for i in range(0, self.meshBCID - len(self.FDSTempFields) - 1):
                self.FDSTempFields.append(None)
            self.FDSTempFields.append(loadedField)
        else:
            self.FDSTempFields[self.meshBCID - 1] = loadedField
        return loadedField

    def loadDumpedInfo(self, filename, filenames):
        text_file = open(filename, "r")
        nFiles = int(text_file.readline())
        for i in range(0, nFiles):
            filenames.append(text_file.readline())
            filenames[i] = filenames[i][:-1]
        text_file.close()
        return nFiles

    def dumpASTField(self, filename="FDSASTF"):
        if self.nVTU_last == self.nVTU:
            self.nVTU = self.nVTU + 1
        text_file = open("%s_s%d.pf" % (filename, self.nVTU), "w")
        text_file.write("%f\n%d\n" % (self.FDSASTField.time, len(self.FDSASTField.value)))
        for i in range(0, len(self.FDSASTField.value)):
            text_file.write("%f\n" % self.FDSASTField.value[i])
        text_file.close()
        text_file = open("%s.tpf" % (filename), "w")
        text_file.write("%d\n" % self.nVTU)
        for i in range(0, self.nVTU):
            text_file.write("%s_s%d.pf\n" % (filename, i + 1))
        text_file.close()

    def loadASTField(self, filename):
        if self.ASTMesh == None:
            filenameM = "./ASTPointsCoords.txt"
            if os.path.exists(filenameM):
                nodeid = 0
                print("loading AST mesh points...")
                newmesh = Mesh.UnstructuredMesh()
                text_file = open(filenameM, "r")
                textLine = text_file.readline()
                while (len(textLine) > 1):
                    splitted = textLine.split()
                    if len(splitted) == 3:
                        nodeid = nodeid + 1
                        newmesh.vertexList.append(
                            Vertex.Vertex(nodeid, nodeid, [float(splitted[0]), float(splitted[1]), float(splitted[2])]))
                    textLine = text_file.readline()
                text_file.close()
                self.ASTMesh = newmesh
                print("%d ASTPoints" % nodeid)

        if not self.ASTMesh:
            print("ASTMesh was not loaded!!!\n")
        # print(self.ASTMesh.vertexList)
        text_file = open(filename, "r")
        stepTime = float(text_file.readline())
        numberOfNodes = int(text_file.readline())

        tempValues = []
        for i in range(0, numberOfNodes):
            tempValues.append(float(text_file.readline()))

        self.FDSASTField = Field.Field(self.ASTMesh, FieldID.FID_Temperature, ValueType.ValueType.Scalar, 'C', stepTime,
                                       tempValues)
        return self.FDSASTField

    def giveTempFieldValues(self):
        return self.FDSTempFields[self.meshBCID - 1].value

    def giveASTFieldValues(self):
        return self.FDSASTField.value

    def giveASTG(self, ind):
        loc_ind = c_int(ind)
        ASTG = c_double(0)
        self.libfds.give_me_astg(byref(loc_ind), byref(ASTG))
        return ASTG.value

    def printASTG(self, ind):
        print(self.giveASTG(ind))


# ##########################################################################################################
# #####    ####    ###     ##     ##  ####  #######   ##   ##   ############################################
# ####  ##  ##  ##  ##  #####  #####   ##   #######   ##   ##   ############################################
# ####  ##  ##  ##  ##     ##    ###        #######  ###  ###  #############################################
# ####  ##  ##  ##  ##  #####  #####  #  #  ################################################################
# #####    ####    ###  #####     ##  ####  ######  ###  ###  ##############################################
# ##########################################################################################################
class oofem_api(Application.Application):
    def __init__(self, inputfilename):
        Application.Application.__init__(self)
        global oofemImported
        global liboofem
        if not oofemImported:
            oofemImported = True
            import liboofem as liboofem
        self.TFBC = TempField()

        self.ASTField = TempField()

        print("OOFEM reading input file '%s'" % (inputfilename))
        super(oofem_api, self).__init__(inputfilename, "./")
        self.dt = 0.0
        self.lastT = 0.0
        self.filename_oofem = inputfilename
        dr = liboofem.OOFEMTXTDataReader(self.filename_oofem)
        self.transpModel = liboofem.InstanciateProblem(dr, liboofem.problemMode._processor, 0)
        self.transpModel.checkProblemConsistency()
        print(self.transpModel.giveClassName())
        activeMStep = self.transpModel.giveMetaStep(1);
        self.transpModel.initMetaStepAttributes(activeMStep)
        self.fieldMan = self.transpModel.giveContext().giveFieldManager()

    def solveStep(self, acttime):
        self.lastT = acttime
        self.transpModel.preInitializeNextStep()
        self.transpModel.giveNextStep()
        previousTargetTime = self.transpModel.givePreviousStep().giveTargetTime()
        currentStep = self.transpModel.giveCurrentStep()
        deltaToofem = acttime - previousTargetTime
        alpha = (currentStep.giveIntrinsicTime() - previousTargetTime) / (
                    currentStep.giveTargetTime() - previousTargetTime)
        currentStep.setTargetTime(previousTargetTime + deltaToofem)
        currentStep.setTimeIncrement(deltaToofem)
        currentStep.setIntrinsicTime(previousTargetTime + alpha * deltaToofem)
        self.transpModel.initializeYourself(currentStep)
        print("OOFEM TimeStep", currentStep.giveNumber(), "target time", currentStep.giveTargetTime(), "alpha", alpha)
        self.transpModel.solveYourselfAt(currentStep)
        self.transpModel.updateYourself(currentStep)
        self.transpModel.terminate(currentStep)

    def complete(self):
        timeStep = self.transpModel.giveCurrentStep()

    def regTempField(self):
        self.fieldMan.registerField(self.TFBC.TField, liboofem.FieldType.FT_TemperatureAmbient)

    def regASTField(self):
        self.fieldMan.registerField(self.ASTField.TField, liboofem.FieldType.FT_TemperatureAmbient)

    def setTempField(self, temp_field_set):
        self.TFBC.setValues(temp_field_set)

    def setASTField(self, temp_field_set):
        self.ASTField.setValues(temp_field_set)

    def createUniformGridTempBCField(self, lo, hi, div, temp=0.):
        self.TFBC.setType(1)
        self.TFBC.TField.setGeometry(lo, hi, div)
        self.TFBC.nnodes = (div[0] + 1) * (div[1] + 1) * (div[2] + 1)
        self.TFBC.ncells = (div[0]) * (div[1]) * (div[2])
        oField = liboofem.FloatArray(self.TFBC.nnodes)
        for i in range(0, self.TFBC.nnodes):
            oField[i] = float(temp)
        self.TFBC.TField.setValues(oField)
        self.fieldMan.registerField(self.TFBC.TField, liboofem.FieldType.FT_TemperatureAmbient)