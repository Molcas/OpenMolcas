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
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: Coord(3,nAtoms)
real(kind=wp), intent(out) :: W1(3,nAtoms*8)
integer(kind=iwp), intent(out) :: iAll_Atom
#include "Molcas.fh"
integer(kind=iwp) :: iAtom, iChAtom, iCo, iCoSet(0:7,0:7), iGen(3), iStab(0:7), MaxDCR, nCoSet, nGen, nStab
integer(kind=iwp), external :: iChxyz

!----------------------------------------------------------------------*
W1(:,1:nAtoms) = Coord(:,:)
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
    call OA(iCoSet(iCo,0),W1(:,iAtom),W1(:,iAll_Atom))

  end do

end do

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Expand_Coor
