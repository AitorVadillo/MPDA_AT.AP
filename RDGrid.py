"""A simulation of two virtual chemicals reacting and diffusing on a 2D grid using the Gray-Scott model.
    Inputs:
        Width: Int, Number of cells in width
        Heigth: Int, Number of cells in Heigth
        Run: Boolean, True to initialize the reaction
    Output:
        a: The a output variable"""
# http://karlsims.com/rd.html
"""
Reaction: two Bs convert an A into B, as if B reproduces using A as food.
Some typical values used, for those interested, are: DA=1.0, DB=.5, f=.055, k=.062 (f and k vary for different patterns), and Δt=1.0. 
The Laplacian is performed with a 3x3 convolution with center weight -1, adjacent neighbors .2, and diagonals .05. 
The grid is initialized with A=1, B=0, and a small area is seeded with B=1.
The grid is visualized by assigning each cell a color from it's A and B values. Here A is white and B is black.
"""
import Rhino.Geometry as rg
import random as rnd
import ghpythonlib.treehelpers as th


# http://karlsims.com/rd.html
# Reaction difusion Variables
a0 = 1        # Initial concentration of substance A In index 0 in the con list
b0 = 0        # Initial concentration of substance B In index 1 in the con list
da = 1       # Diffusion rate for substance A
db = 0.5     # Diffusion rate for substance B
f = 0.055    # Feed rate Speed of adding substance A
k = 0.062    # Kill rate Speed of removing substance B
dt = 1       # Delta time for each iteration


mat0 = []    # Matrix with the concentration stets in the previous iteration or initial state (List of Columns)
mat = []     # Matrix with the concentration states (List of Columns)
col = []     # List of column
temp = []    # Temporal list

# Creates the swap between matrix mat goes to ma0 in order to calculate mat
def changematrix():
    #temp = mat0
    mat0 = mat
    #mat = temp
    
# Calculates the Laplacian Convolution for substance  A
def laplacianA ():
    resA = 1
    return resA

# Calculates the Laplacian Convolution for substance B
def laplacianB ():
    resB = 1
    return resB

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
# Creates the concentration in mat0
for x in range(3,6):
    for y in range(3,6):
        mat[x][y][1] = 1
            

##################################################

# Iteration loop
while Run == True:
    
    # Calculates the new state
    for x in range(1, Width-1):
        for y in range(1, Heigth-1):
            a0 = mat0[x][y][0]
            b0 = mat0[x][x][1]
            
            a = a0 + (da * laplacianA() - a0*(b0**2) + f * (1 - a0) ) * dt
            b = b0 + (db * laplacianB() + a0*(b0**2) - (k + f) * b0 ) * dt
            
            mat[x][y][0] = limiter(a, 0, 1)
            mat[x][y][1] = limiter(b, 0, 1)
           
            
    
    
            
    Run = False





######DEBUGING######
c = th.list_to_tree(mat)
d = th.list_to_tree(mat0)
mesh = rg.Mesh.CreateFromPlane(rg.Plane(rg.Point3d(0,0,0), rg.Vector3d.ZAxis), xSize, ySize, Width, Heigth)
print("Finish")


