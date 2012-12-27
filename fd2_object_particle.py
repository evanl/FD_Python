import numpy as np


class ParticleField(object):

  def __init__(xp_in, yp_in, headMatrix, xLocations, yLocations,\
               dx, dy, nx, ny, conductivity, porosity):
    """
    xp_in and yp_in are initial coordinates of particle
    headMatrix is a numpy array of head values
    where xLocations and yLocations correspond to the grid centers
    dx and dy are grid spacing.
    nx and ny are stored for convenience
    conductivity and porosity are constant across the domain in this version.
    """
    self.xp = [xp_in]
    sely.yp = [yp_in]
    self.head = headMatrix
    self.xgrid = xLocations
    self.ygrid = yLocations
    self.dx = dx
    self.dy = dy
    self.nx = nx
    self.ny = ny
    self.k = conductivity
    self.phi = porosity
    self.t = [0.]


  def ParticleTrack(tMax):
    """ 
    Runs the particle tracking for the parameters that were initialized. 
    Follows the pollock algorithm prett closely this time.
    """
    outOfDomain = False

    # Assign particle to cell 
    # assigns the index of the cell being used in the arrays. 
    AssignInitialCell()

    while outOfDomain == False and self.t < tMax:

      # compute cell face velocities
      GridCell.CalcFaceVelocities(self.head[i,j], self.head[i-1,j]\
          self.head[i+1,j], self.head[i,j-1], self.head[i,j+1])

      # determine exit faces

      # Compute transit time, determine actual exit face

      # compute exit point coordinates

      # write particle coordinates and times to list
      self.xp.append(xp_out)
      self.yp.append(yp_out)

      # Determine new cell
      AssignNewCell(exitFace )

class GridCell(xCenter, yCenter, dx, dy, conductivity, porosity):
  """ Grid cell contains the edge coordinates, face velocities and 
  if the particle is a capture point. 
  Takes in the center and cell widths in the following form:
  xCenter, yCenter, dx, dy
  and stores the edge faces as a 2d list.
  x0[0] = xwest
  x0[1] = ysouth
  x1[0] = xeast
  x1[1] = ynorth
  Also stores the porosity and conductivity for transmissivity equations.
  """
  self.x0 = [xCenter - dx/2, yCenter - dy/2]
  self.x1 = [xCenter + dx/2, yCenter + dy/2]
  self.exit0 = [True, True]
  self.exit1 = [True, True]
  self.dx = [dx, dy]
  self.v0 = [0., 0.]
  self.v1 = [0., 0.]
  self.A = [0., 0.]
  self.k = conductivity
  self.phi = porosity

  def CalcFaceVelocities(hC, hW, hE, hS, hN):
    """ 
    Calculates the velocities on each face by linearly interpolating between
    the heads of the 5-point finite difference
    """
    
    self.v0[0] = k * (hW - hC) / (self.dx[0] * self.phi)
    self.v0[1] = k * (hC - hE) / (self.dx[1] * self.phi)
    self.v1[0] = k * (hS - hC) / (self.dx[0] * self.phi)
    self.v1[1] = k * (hC - hN) / (self.dx[1] * self.phi)
    
    for i in range(2):
      A[i] = (v1[i] - v0[i])/ dx[i]
    
  def AssignExitFaces(

   
  
