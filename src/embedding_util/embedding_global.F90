!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
! Data shared for an embedding potential on a grid

module Embedding_Global

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: nEmbGridPoints
! The energy coming from the embedding potential
real(kind=wp) :: Eemb

!****** Arrays
! Coordinates of the grid points, stored as x1, y1, z1, x2, y2...
real(kind=wp), allocatable :: embGridCoord(:,:)
! The actual value of the embedding potential at grid point i
real(kind=wp), allocatable :: embPotVal(:)
! The corresponding weighting factor of grid point i
real(kind=wp), allocatable :: embWeight(:)
! The integrals in form <a|embPot|b>
real(kind=wp), allocatable :: embInt(:)

!****** Flags
! Flag whether an embedding potential is used
logical(kind=iwp) :: embPot
! Flag whether the potential is given in a basis representation
logical(kind=iwp) :: embPotInBasis
! Flag whether a path to an output grid has been specified,
! i.e. whether the output grid is different from the input grid
logical(kind=iwp) :: outGridPathGiven
! What output is requested
logical(kind=iwp) :: embWriteDens, embWriteEsp, embWriteGrad, embWriteHess
logical(kind=iwp) :: embDebug

!****** File paths
! The path to the file containing grid and embedding potential
character(len=256) :: embPotPath
! The output grid (if different from input)
character(len=256) :: outGridPath
! The output data, for density, electrostatic potential and its gradient
character(len=256) :: embOutDensPath, embOutEspPath, embOutGradPath, embOutHessPath

public :: nEmbGridPoints, Eemb, embGridCoord, embPotVal, embWeight, embInt, embPot, embPotInBasis, outGridPathGiven, embWriteDens, &
          embWriteEsp, embWriteGrad, embWriteHess, embDebug, embPotPath, outGridPath, embOutDensPath, embOutEspPath, &
          embOutGradPath, embOutHessPath

end module Embedding_Global
