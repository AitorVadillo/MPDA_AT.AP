"""A simulation of two virtual chemicals reacting and diffusing on a 2D grid using the Gray-Scott model.
    Inputs:
        Width: Int, Number of cells in width
        Heigth: Int, Number of cells in Heigth
        Iter: Int, Number of iterations
        xSize: Int, Distance between points in X direction
        ySize: Int, Distance between points in Y direction
        amp: Int, Amplitude of the curves
        CompuCurves: Bool, If true computs and outputs the curves
        CompuSrf: Bool, If true computs and outputs the curves and the surface 
    Output:
        srf: Srf, The surface
        mesh: Mesh, The mesh
        """
        
# http://karlsims.com/rd.html
"""
A simulation of two virtual chemicals reacting and diffusing on a 2D grid using the Gray-Scott model.
Reaction: two Bs convert an A into B, as if B reproduces using A as food.
Some typical values used, for those interested, are: DA=1.0, DB=.5, f=.055, k=.062 (f and k vary for different patterns), and Δt=1.0. 
The Laplacian is performed with a 3x3 convolution with center weight -1, adjacent neighbors .2, and diagonals .05. 
The grid is initialized with A=1, B=0, and a small area is seeded with B=1.
The grid is visualized by assigning each point a Z coordinate from it's A and B values. Here A is Z=amp  and B Z=-amp.
"""
import Rhino.Geometry as rg
import random as rnd
import ghpythonlib.treehelpers as th


# http://karlsims.com/rd.html
# Reaction difusion Variables
a0 = 1       # Initial concentration of substance A In index 0 in the con list
b0 = 0       # Initial concentration of substance B In index 1 in the con list
da = 1       # Diffusion rate for substance A
db = 0.5     # Diffusion rate for substance B
f = 0.055    # Feed rate Speed of adding substance A
k = 0.062    # Kill rate Speed of removing substance B
dt = 1       # Delta time for each iteration

##################################################
# Creates the matrix of points
Pts = []
colPts = []

for x in range(0, Width, xSize):
    colPts = []
    for y in range(0, Heigth, ySize):
        pt = rg.Point3d(x, y, 0)
        colPts.append(pt)
    Pts.append(colPts)



mat0 = []    # Matrix with the concentration stets in the previous iteration or initial state (List of Columns)
mat = []     # Matrix with the concentration states (List of Columns)
col = []     # List of colum
temp = []    # Temporal list
corZ = []    # Matrix with cordinate Z of Pts
curves = []  # List with interpolated curves

# Creates the swap between matrix mat goes to ma0 in order to calculate mat
def changematrix():
    global mat0 
    global mat
    #temp = mat0
    mat0 = mat 
    #mat = temp
    
# Calculates the Laplacian Convolution for substance  A
def laplacianA (x, y):
    sumA = 0
    sumA += mat0[x-1][y+1][0] * 0.05
    sumA += mat0[x][y+1][0] * 0.2
    sumA += mat0[x+1][y+1][0] * 0.05
    sumA += mat0[x-1][y][0] * 0.2
    sumA += mat0[x][y][0] * -1
    sumA += mat0[x+1][y][0] * 0.2
    sumA += mat0[x-1][y-1][0] * 0.05
    sumA += mat0[x][y-1][0] * 0.2
    sumA += mat0[x+1][y-1][0] * 0.05
    return sumA

# Calculates the Laplacian Convolution for substance B
def laplacianB (x, y):
    sumB = 0
    sumB += mat0[x-1][y+1][1] * 0.05
    sumB += mat0[x][y+1][1] * 0.2
    sumB += mat0[x+1][y+1][1] * 0.05
    sumB += mat0[x-1][y][1] * 0.2
    sumB += mat0[x][y][1] * -1
    sumB += mat0[x+1][y][1] * 0.2
    sumB += mat0[x-1][y-1][1] * 0.05
    sumB += mat0[x][y-1][1] * 0.2
    sumB += mat0[x+1][y-1][1] * 0.05
    return sumB
    
# Returns the number if is between min/max, if not gives the min or the max
def limiter(val, minVal, maxVal):
    return min(maxVal, max(minVal, val))
    
##################################################
# Creates the matrix mat0 and mat
for x in range(Width):
    for y in range(Heigth):
        col.append([a0, b0])
    mat0.append(col)
    col = []

for x in range(Width):
    for y in range(Heigth):
        col.append([a0, b0])
    mat.append(col)
    col = []

##################################################
# Creates the initial concentration in mat0
for x in range(45,56):
    for y in range(45,56):
        mat0[x][y][1] = 1

# Set a=1 and b=1 to the boundary points, in order to layout in z = 0
for x in range(Width):
    for y in range(Heigth):
        if x == 0 or x == Width-1 or y == 0 or y == Heigth-1 :
            mat0[x][y][1] = 1

##################################################
# Set the preconditions in Z (Initial case)
#--------------------------------------------------------------------------------------------------------
# Matrix of points Z Values
for x in range(Width):
    col = []
    for y in range(Heigth):
        z = (mat0[x][y][0] * amp) - (mat0[x][y][1] * amp)
        col.append(z)

    corZ.append(col)
    
    # Set the new value for Z on the Pts
for x in range(Width):
    for y in range(Heigth):
        Pts[x][y].Z = corZ[x][y]
"""
# List of interpolated curves
if CompuSrf == True or CompuCurves == True:
    for i in range(len(Pts)):
            curves.append(rg.Curve.CreateInterpolatedCurve(Pts[i], 3))

# Create the surface
if CompuSrf == True:
    srf = rg.Brep.CreateFromLoft(curves, rg.Point3d.Unset, rg.Point3d.Unset, rg.LoftType(), False)

# Create the mesh
meshVertex = []
meshVertexL = []
mesh = rg.Mesh()

for i in range(Width):
    for j in range(Heigth):
        mesh.Vertices.Add(Pts[i][j])

# Function that converts the position in the matrix (r,c) to an index
def rc_to_index (co, ro):
    index = co * Heigth + ro
    return index

for i in range(1, Width):
        
    for j in range(0, Heigth-1):
        mesh.Faces.AddFace(rc_to_index(i, j), rc_to_index(i, j+1), rc_to_index(i-1, j+1), rc_to_index(i-1, j))
"""
#--------------------------------------------------------------------------------------------------------

##################################################
# Iteration loop
for i in range(Iter):
    curves = []
    corZ = []
    # Calculates the new state
    for x in range(1, Width-1):
        for y in range(1, Heigth-1):
            a0 = mat0[x][y][0]
            b0 = mat0[x][x][1]
            
            a = a0 + (da * laplacianA(x, y) - a0* (b0**2) + f * (1 - a0) ) * dt
            b = b0 + (db * laplacianB(x, y) + a0* (b0**2) - (k + f) * b0 ) * dt
            
            mat[x][y][0] = limiter(a, 0, 1)
            mat[x][y][1] = limiter(b, 0, 1)
            
    changematrix()

##################################################
#Calculate geometry after iterations
#Remove Case 0
 


# Matrix of points Z Values
for x in range(Width):
    col = []
    for y in range(Heigth):
            z = (mat[x][y][0] * amp) - (mat[x][y][1] * amp)
            col.append(z)
    corZ.append(col)
    
# Set the new value for Z on the Pts
for x in range(Width):
    for y in range(Heigth):
        Pts[x][y].Z = corZ[x][y]
    
# Create the nurbs per column
if CompuSrf == True or CompuCurves == True:
    for i in range(len(Pts)):
        curves.append(rg.Curve.CreateInterpolatedCurve(Pts[i], 3))
    
# Create the surface
if CompuSrf == True:
    srf = rg.Brep.CreateFromLoft(curves, rg.Point3d.Unset, rg.Point3d.Unset, rg.LoftType(), False)
    
# Create the mesh
meshVertex = []
meshVertexL = []
mesh = rg.Mesh()

for i in range(Width):
    for j in range(Heigth):
        mesh.Vertices.Add(Pts[i][j])

# Function that converts the position in the matrix (r,c) to an index
def rc_to_index (co, ro):
    index = co * Heigth + ro
    return index

for i in range(1, Width):
        
    for j in range(0, Heigth-1):
        mesh.Faces.AddFace(rc_to_index(i, j), rc_to_index(i, j+1), rc_to_index(i-1, j+1), rc_to_index(i-1, j))

######DEBUGING######
#previous = th.list_to_tree(mat0)
#final = th.list_to_tree(mat)
Grid = th.list_to_tree(Pts)
print("Finish")
print (mat0)