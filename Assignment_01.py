import Rhino.Geometry as rg
import math
import ghpythonlib.treehelpers as th

"""Provides a scripting component.
    Inputs:
        x: Int, Extension in x direction of the mesh
        y: Int, Extension in y direction of the mesh
    Output:
        mesh: Mesh, The resulting ondulated mesh"""


# Points lists
ptsUp = []
ptsDw = []

# Filling the points lists with points
for i in range(x):
    ptsDw.append(rg.Point3d(i, 0, 0))
    ptsUp.append(rg.Point3d(i, y, 0))
# Debugging output
a = ptsUp
b = ptsDw
#---------------------------------------------------------------

# Lines list
lines = []


# Filling the lines list with lines reparametericed
for i in range(len(ptsUp)):
    crv = (rg.LineCurve(ptsDw[i], ptsUp[i]))
    crv.Domain = rg.Interval(0, 1)
    lines.append(crv)
    
# Debugging output
c = lines
#---------------------------------------------------------------

# Curve points list of lists Pts contains linePts lists
Pts = []
linePts = []

for i in range(len(lines)):
    linePts = []
    for j in range(0, y+1):
        linePts.append(lines[i].PointAt(j/y))
    Pts.append(linePts)

# Debugging output
d = th.list_to_tree(Pts)
#---------------------------------------------------------------

# Zets list of lists zets contains zetsL lists (Mantain strucutre Why not :D)
zets = []
zetsL = []

for i in range(len(Pts)):
    zetsL = []
    for j in range(len(Pts[i])):
        vec = rg.Vector3d(Pts[i][j].X, Pts[i][j].Y, 0)
        Len = vec.Length
        zet = math.sin(Len)
        zetsL.append(zet)
    zets.append(zetsL)

# Debugging output
e = th.list_to_tree(zets)
#---------------------------------------------------------------

# Copy Pts matrix to move the points
PtsC = []
PtsLC = []

for i in range(len(Pts)):
    PtsLC = []
    for j in range(len(Pts[i])):
        PtsLC.append(Pts[i][j])
    PtsC.append(PtsLC)

# Moved points list moved of lists contains movedL lists
moved = []
movedL = []

for i in range(len(PtsC)):
    movedL = []
    for j in range(len(PtsC[i])):
        PtsC[i][j].Z = zets[i][j]
        movedL.append(PtsC[i][j])
    moved.append(movedL)
    
# Debugging output
f = th.list_to_tree(moved)
#---------------------------------------------------------------

# List with interpolated crvs
crvs = []
for i in range(len(moved)):
    crvs.append(rg.Curve.CreateInterpolatedCurve(moved[i], 3))

# Debugging output  
g = crvs
#---------------------------------------------------------------

# Creating the surface
srf = rg.Brep.CreateFromLoft(crvs, rg.Point3d.Unset, rg.Point3d.Unset, rg.LoftType(), False)

# Debugging output  
h = srf
#---------------------------------------------------------------

# Creating the mesh
meshVertex = []
meshVertexL = []
mesh = rg.Mesh()

for i in range(len(moved)):
    for j in range(len(moved[i])):
        mesh.Vertices.Add(moved[i][j])

# Function that converts the position in the matrix (r,c) to an index
def rc_to_index (co, ro):
    index = co * (len(moved[0])) + ro
    return index

for i in range(1, len(moved)):
    
    for j in range(0, len(moved[i])-1):
        mesh.Faces.AddFace(rc_to_index(i, j), rc_to_index(i, j+1), rc_to_index(i-1, j+1), rc_to_index(i-1, j))


#mesh.Normals.ComputeNormals()
#mesh.Compact()




#mesh = rg.Mesh.CreateFromBrep(srf[0], rg.MeshingParameters(1,1))
#mesh = rg.Mesh.QuadRemeshBrep(srf[0], rg.QuadRemeshParameters())