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

module grid_it_globals

use Definitions, only: wp, iwp

implicit none
private

#include "LenIn.fh"

! GridSparse, GridDense number of points in a.u.
! MAXGRID maximum number of generated grid ONLY if you set up it by hand!
integer(kind=iwp), parameter :: GridDense = 10, GridNormal = 3, GridSparse = 2, MAXGRID = 100
character(len=*), parameter :: VERSION = '1.01'

! IMaxup/Down 0 default number of orbitals near HOMO-LUMO
! iGridNpt - 3 integer to set GRID
!            in Natural way i.e. '10' means at 0.0 0.1 0.2 etc.
! iBinary  - create binary file, not ASCII to output GRID
! TheGap   - distance between molecule and automatically generated box
integer(kind=iwp) :: iGauss, iGridNpt(3), iMaxDown, iMaxUp, iReq(MAXGRID*2), iAuMO, iBinary, iDerivative, itRange, nAtoms, &
                     nBytesPackedVal, nGridPoints, nReq, LevelPrint, LID, LID_ab, LuVal, LuVal_ab
logical(kind=iwp) :: isAll, isAtom, isColor, isCurDens, isCutOff, isDebug, isDensity, isLine, isLuscus, isMOPack, isSphere, &
                     isTheOne, isTotal, isUHF, isUserGrid, isVirt, NoOrb, NoSort
real(kind=wp) :: CutOff, GridAxis1(3), GridAxis2(3), GridAxis3(3), GridOrigin(3), OneCoor(7), Region(2), TheGap, Virt
real(kind=wp), allocatable :: Coor(:,:), Grid(:,:)
character(len=LenIn), allocatable :: AtomLbl(:)
character(len=256) :: TheName
character(len=80) :: Title1

public :: AtomLbl, Coor, CutOff, Grid, GridAxis1, GridAxis2, GridAxis3, GridDense, GridNormal, GridSparse, GridOrigin, iGauss, &
          iGridNpt, iMaxDown, iMaxUp, iReq, isAll, isAtom, iAuMO, iBinary, isColor, isCurDens, isCutOff, isDebug, isDensity, &
          iDerivative, isLine, isLuscus, isMOPack, isSphere, isTheOne, isTotal, isUHF, isUserGrid, isVirt, itRange, LevelPrint, &
          LID, LID_ab, LuVal, LuVal_ab, MAXGRID, nAtoms, nBytesPackedVal, nGridPoints, NoOrb, NoSort, nReq, OneCoor, Region, &
          TheGap, TheName, Title1, VERSION, Virt

end module grid_it_globals
