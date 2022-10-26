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

subroutine Expand_Coor(Coord,nAtoms,W1,iAll_Atom)
!***********************************************************************
!                                                                      *
!     purpose: to generate a list of all atoms                         *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: iOper, nIrrep

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
real*8 Coord(3,nAtoms)
real*8 W1(3,nAtoms*8)
integer iGen(3), iCoSet(0:7,0:7), iStab(0:7)

!----------------------------------------------------------------------*
call dcopy_(nAtoms*3,Coord,1,W1,1)
!----------------------------------------------------------------------*
! Apply the symmetry operations                                        *
!----------------------------------------------------------------------*
nGen = 0
if (nIrrep == 2) nGen = 1
if (nIrrep == 4) nGen = 2
if (nIrrep == 8) nGen = 3
if (nGen >= 1) iGen(1) = iOper(1)
if (nGen >= 2) iGen(2) = iOper(2)
if (nGen >= 3) iGen(3) = iOper(4)

MaxDCR = 0
iAll_Atom = nAtoms
do iAtom=1,nAtoms
  iChAtom = iChxyz(W1(1,iAtom),iGen,nGen)
  call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
  nCoSet = nIrrep/nStab

  do iCo=1,nCoSet-1

    iAll_Atom = iAll_Atom+1
    call OA(iCoSet(iCo,0),W1(1:3,iAtom),W1(1:3,iAll_Atom))

  end do

end do

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Expand_Coor
