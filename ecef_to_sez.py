# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#  Converts ECEF vector components to SEZ 
#  See "Fundamentals of Astrodynamics and Applications, Fourth Edition" by
#  David A. Vallado, pages 172-173
# Parameters:
#  o_x_km: ECEF origin x-component in km
#  o_y_km: ECEF origin y-component in km
#  o_z_km: ECEF origin z-component in km
#  x_km: ECEF x-component in km
#  y_km: ECEF y-component in km
#  z_km: ECEF z-component in km
# Output:
#  Prints the sez components from the ECEF Origin
#
# Written by Tejas Vinod
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.
#
# Test Inputs:
# python3 ecef_to_sez.py 822.933 -4787.187 4120.262 1131.698 -4479.324 4430.228
# -398.778, 356.458, 10.340
# 
# 

# import Python modules
import math # math module
import sys  # argv

# "constants"
R_E_KM = 6378.137
E_E    = 0.081819221456

# helper functions

## calculated denominator
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-(ecc**2)*(math.sin(lat_rad)**2))

## matrix multiplication
def matrix_multiply(matrix1, matrix2):
    # Get the dimensions of the matrices
    rows_matrix1 = len(matrix1)
    cols_matrix1 = len(matrix1[0])
    rows_matrix2 = len(matrix2)
    cols_matrix2 = len(matrix2[0])
    
    # Check if matrices are able to be multiplied
    if cols_matrix1 != rows_matrix2:
        raise ValueError("Matrix multiplication is not possible. Columns of matrix1 must equal rows of matrix2.") 
    
    # Perform the multiplication
    result = [[0 for _ in range(cols_matrix2)] for _ in range(rows_matrix1)]
    for i in range(rows_matrix1):
        for j in range(cols_matrix2):
            for k in range(cols_matrix1):
                result[i][j] += matrix1[i][k] * matrix2[k][j]
    
    return result

# initialize script arguments
o_x_km = float('nan') # ECEF origin x-component in km
o_y_km = float('nan') # ECEF origin y-component in km
o_z_km = float('nan') # ECEF origin z-component in km
x_km = float('nan') # ECEF x-component in km
y_km = float('nan') # ECEF y-component in km
z_km = float('nan') # ECEF z-component in km

# parse script arguments
if len(sys.argv)==7:
  o_x_km = float(sys.argv[1])
  o_y_km = float(sys.argv[2])
  o_z_km = float(sys.argv[3])
  x_km = float(sys.argv[4])
  y_km = float(sys.argv[5])
  z_km = float(sys.argv[6])
else:
  print(\
   'Usage: '\
   'python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km'\
  )
  exit()

# write script below this line

# Find Lat Lon of SEZ Origin
# calculate longitude
o_lon_rad = math.atan2(o_y_km,o_x_km)

# initialize lat_rad, r_lon_km, r_z_km
o_lat_rad = math.asin(o_z_km/math.sqrt(o_x_km**2+o_y_km**2+o_z_km**2))
r_lon_km = math.sqrt(o_x_km**2+o_y_km**2)
prev_lat_rad = float('nan')

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(o_lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,o_lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = o_lat_rad
  o_lat_rad = math.atan((o_z_km+c_E*(E_E**2)*math.sin(o_lat_rad))/r_lon_km)
  count = count+1

# Find SEZ
r_rel_ecef = [[x_km-o_x_km], [y_km-o_y_km], [z_km-o_z_km]]

Ryinv = [ 
      [math.sin(o_lat_rad), 0, -math.cos(o_lat_rad)],
      [                   0, 1,                   0],
      [math.cos(o_lat_rad), 0,  math.sin(o_lat_rad)] ] # ry(90-phi)^-1
Rzinv = [ 
      [ math.cos(o_lon_rad), math.sin(o_lon_rad), 0],
      [-math.sin(o_lon_rad), math.cos(o_lon_rad), 0],
      [                  0,                    0, 1] ] # rz(theta)^-1

r_sez = matrix_multiply(Ryinv, matrix_multiply(Rzinv, r_rel_ecef))
    
print(r_sez[0][0])
print(r_sez[1][0])
print(r_sez[2][0])   