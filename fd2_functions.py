from time import time, clock
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as splinalg

class BoundaryCondition(object):
    
  def __init__(self, w_t = "head", w_v = 0,  \
                       e_t = "head", e_v = 0, \
                       s_t = "head", s_v = 0, \
                       n_t = "head", n_v = 0 ):

        # if a boundary condition value is entered as a function, make sure to change boundary condition routine in SpatialSolve
    self.west_type = w_t
    self.west_val = w_v

    self.east_type = e_t
    self.east_val = e_v

    self.south_type = s_t
    self.south_val = s_v

    self.north_type = n_t
    self.north_val = n_v
#########################
def Transmissivity(x, y, B = 10.0, xcond = 0.001, ycond = 0.001):
  def thickness(x, y):
    x = x * 3.281
    y = y * 3.281 
    if x < 550.:
      bottom = 30. - (30 - -50.)/(850. - 0) * x
    elif x < 850.:
      bottom = -50
    elif x < 1000.:
      bottom = -50. + (10 - -50)/(1000.-850) * (x-850.)
    elif x < 2250.:
      bottom = 10 + (70. - 10) / (2250. - 1000) * ( x-1000.)
    else :
      bottom = 70.

    if x < 1000:
      top = 50
    else :
      top = 50 + (90 - 50) / (2250. - 1000) * (x - 1000.)

    # returns it in meters
    return (top - bottom)/3.281

  # need conductivity function
  def conductivity(x,y):
    x = x * 3.281
    if x < 1350:
      k = 100.
    else : 
      k = 10.
    return k / 3.281
  
  B = thickness(x,y)

  xcond = conductivity(x,y)
  ycond = xcond
  
  T = []

  T.append(B*xcond)
  T.append(B*ycond)

  return T
##########################
def Source(x,y,t, dx = 1., dy = 1., \
    infil = False, infilval = 0., \
    well = False, Qwell = [0.0], wellLoc = [[0,0]] ): 
    
  source = 0.0
    #infiltration cases
  if infil != False:
    source = infilval

  # well case 
  if well != False:
    for i in range(len(wellLoc)):
      if (wellLoc[i][0] - dx/2)<  x  < (wellLoc[i][0] + dx/2):
        if  (wellLoc[i][1] - dy/2)<  y  < (wellLoc[i][1] + dy/2):
          print "WEll " + str(i) + " IS ON"
          source = source + Qwell[i] / (dx * dy)

  return source
##########################
def HarmAvg(T1, T2, dx1, dx2):
    return (dx1 + dx2) / ( (dx1/T1) + (dx2/T2))
##########################
def SpatialSolve(xw, xe, ys, yn, nx, ny, dx, dy, bc, \
    t = 0, S = 0, dt = 1.0, r= np.zeros((1,1)), h0 = 0.0, K_B =0.0,\
    wells=False, Qwell1=[0.0], wellLoc1=[[0,0]],
    infil = False, infilval = 0): 
    
  # # # 
  # Spatial Solve 
  # Returns a 2d array of average head values for vertically 
  # integrated groundwater flow. 
  # 
  # OUTPUT
  N = (nx-2) * (ny -2)

  Loff = np.zeros((N,1))
  L    = np.zeros((N,1))
  R    = np.zeros((N,1))
  Roff = np.zeros((N,1))

  if S == 0: 
    #A = np.zeros((N, N))
    C = np.zeros((N,1))
  else:
    #A =  S / dt * np.eye(N) 
    C = S/dt * np.ones((N,1))

  if r.all() == 0:    
    r = np.zeros((N,1))

  tinit = time()

  for i in xrange(N):
    
    # operations to determine the center node's global x and y coordinates. 
    col = ( i % (nx-2) ) + 1        
    row = i // (nx-2) + 1 
    if col == 0 :
      col = (nx-2)

    x = xw + dx * col
    y = ys + dy * row

    r[i] += S / dt
    r[i] += dx * dy * Source(x,y,t,dx, dy, well = wells,Qwell= Qwell1,wellLoc =  wellLoc1,
        infil= infil, infilval= infilval)

    #leakage
    if h0 != 0.0:
      # river leakage stuff for final project problem
      xriver = 870/3.281
      if (xriver - dx/2)<  x  < (xriver + dx/2):
        h0 = (45-43)/(3.281 * (yn - ys)) * y + 43/3.281
        r[i] -= K_B * h0 * dx * dy
        C[i] -= K_B * dx * dy
      
    TC = Transmissivity(x     , y     )
    TW = Transmissivity(x - dx, y     )
    TE = Transmissivity(x + dx, y     ) 
    TS = Transmissivity(x     , y - dy)
    TN = Transmissivity(x     , y + dy)

    CW = 2.0 * dy / (dx + dx) * HarmAvg(TC[0], TW[0], dx, dx)
    CE = 2.0 * dy / (dx + dx) * HarmAvg(TC[0], TE[0], dx, dx)
    CS = 2.0 * dx / (dy + dy) * HarmAvg(TC[1], TS[1], dy, dy)
    CN = 2.0 * dx / (dy + dy) * HarmAvg(TC[1], TN[1], dy, dy)

    # am i on west face
    if  col == 1 :
      if bc.west_type == "flux":
        r[i] -= bc.west_val
      else:
        r[i] -= CW * bc.west_val
        #A[i,i] -= CW
        C[i] -= CW
    else:
      #A[i, i-1] = CW
      L[i-1] = CW
      #A[i,i] -= CW
      C[i] -= CW

    # east face
    if col == (nx-2):
      if bc.east_type == "flux":
        r[i] += bc.east_val
      else:
        r[i] -= CE * bc.east_val
        #A[i,i] -= CE
        C[i] -= CE
    else:
      #A[i,i+1] = CE
      R[i+1] = CE
      #A[i,i] -= CE
      C[i] -= CE

    # south face
    if i < (nx-2):
      if bc.south_type == "flux":
        r[i] -= bc.south_val
      else:
        r[i] -= CS * bc.south_val
        #A[i,i] -= CS
        C[i] -= CS
    else:
      #A[i, i-(nx-2)] = CS
      Loff[i-(nx-2)] = CS
      #A[i,i] -= CS
      C[i] -= CS

    # north face
    if i >= (N - (nx-2)):
      if bc.north_type == "flux":
        r[i] += bc.north_val
      else:
        r[i] -= CN * bc.north_val
        #A[i,i] -= CN
        C[i] -= CN
    else:
      #A[i,i +(nx-2)] = CN
      Roff[i+(nx-2)] = CN
      #A[i,i] -= CN
      C[i] -= CN

  # end of  loop
  tpop = time()
  print "matrix create time = " +str(tpop - tinit)

  # create new sparse matrix from pieces 
  data = np.hstack((Loff, L, C, R, Roff))
  data = np.transpose(data)
  diags = np.array([-(nx-2), -1, 0, 1, (nx-2)])
  Abuilt = sp.spdiags(data,diags,N,N)

  # convert to sparse matrix for solve
  #Asparse = sp.dia_matrix(A)
  
  #print Abuilt.todense()
  b = splinalg.spsolve(Abuilt,r)

  tmatsolve = time()
  print "linalg solve time=  " + str( tmatsolve - tpop)
  
  
  H = np.zeros((ny,nx))

  # enter solved values into spatial head matrix
  for i in xrange(1,nx-1):
      for j in xrange(1,ny-1):
          H[j,i] = b[ (i-1) + (nx-2) * (j-1) ]

  # note that boundary condition routine must be changed if the values are functions. 

  for i in range(ny):
      if bc.west_type == "head":
          H[i,0] = bc.west_val            
      else:
          H[i,0] = H[i,1] + 1.0/dx * bc.west_val

  for i in range(ny):
      if bc.east_type == "head":
          H[i,-1] = bc.east_val
      else:
          H[i,-1] = H[i,-2] - 1.0/dx * bc.east_val

  for i in xrange(0,nx):
      if bc.south_type == "head":
          H[0,i] = bc.south_val
      else:
          H[0,i] = H[1,i] + 1.0/dy * bc.south_val

  for i in xrange(0,nx):
      if bc.north_type == "head":
          H[-1,i] = bc.north_val
      else:
          H[-1,i] = H[-2,i] + 1/dy * bc.south_val

  return H
###########################
def TimeSolve( xw , xe, ys, yn, nx, ny, dx, dy, bc, tmax=1, nt=1, S=0, h0 = 0, K_B = 0):

    N = (nx-2) * (ny-2)

    r = np.zeros((N,1))

    print tmax
    print nt
    dt = tmax / (nt-1)


    # initial condition
    
    H = np.zeros((ny,nx,nt))
    print np.shape(H)
    H[:,:,0] = SpatialSolve(xw, xe, ys, yn, nx, ny, dx, dy, bc, h0, K_B)

    for k in xrange(1,nt):
        
        count = 0
        for i in xrange(1,nx-1):
            for j in xrange(1,ny-1):
                r[count,0] = S / dt * H[i,j,k-1]
                count += 1
        t = k * dt

        H[:,:,k] = SpatialSolve(xw, xe, ys, yn, nx, ny, dx, dy, bc, t, S, dt, r, h0, K_B)

    return H
