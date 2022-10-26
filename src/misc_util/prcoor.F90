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

subroutine PrCoor()
!***********************************************************************
!                                                                      *
!     purpose: Write coordinates.                                      *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: Symmetry_Info_Get

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "angstr.fh"
integer iGen(3), iCoSet(0:7,0:7), iStab(0:7), iOper(0:7)
character(LEN=LENIN) AtomLbl(MxAtom), Byte4
real*8, allocatable :: W1(:,:)

!----------------------------------------------------------------------*
! Read no.of symm. species                                             *
!----------------------------------------------------------------------*
call Get_iScalar('nSym',nSym)
!----------------------------------------------------------------------*
! Read symm. oper per symm. species                                    *
!----------------------------------------------------------------------*
call Get_iArray('Symmetry operations',iOper,nSym)
!----------------------------------------------------------------------*
! Read no. of unique atoms in the system                               *
!----------------------------------------------------------------------*
call Get_iScalar('Unique atoms',nAtoms)
!----------------------------------------------------------------------*
! Read atom labels                                                     *
!----------------------------------------------------------------------*
call Get_cArray('Unique Atom Names',AtomLbl,LENIN*nAtoms)
!----------------------------------------------------------------------*
! Read coordinates of atoms                                            *
!----------------------------------------------------------------------*
call mma_Allocate(W1,3,8*nAtoms)
call Get_dArray('Unique Coordinates',W1,3*nAtoms)
!----------------------------------------------------------------------*
! Read nuclear repulsion energy                                        *
!----------------------------------------------------------------------*
call Get_dScalar('PotNuc',PotNuc)
!----------------------------------------------------------------------*
! Pick up the symmetry information                                     *
!----------------------------------------------------------------------*
call Symmetry_Info_Get()
!----------------------------------------------------------------------*
! Apply the symmetry operations                                        *
!----------------------------------------------------------------------*
nGen = 0
if (nSym == 2) nGen = 1
if (nSym == 4) nGen = 2
if (nSym == 8) nGen = 3
if (nGen >= 1) iGen(1) = iOper(1)
if (nGen >= 2) iGen(2) = iOper(2)
if (nGen >= 3) iGen(3) = iOper(4)

iAll_Atom = 0
MaxDCR = 0
iAll_Atom = nAtoms
do iAtom=1,nAtoms
  iChAtom = iChxyz(w1(1:3,iAtom),iGen,nGen)
  call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
  nCoSet = nSym/nStab
  Byte4 = AtomLbl(iAtom)

  do iCo=1,nCoSet-1

    iAll_Atom = iAll_Atom+1
    call OA(iCoSet(iCo,0),W1(1:3,iAtom),W1(1:3,iAll_Atom))
    AtomLbl(iAll_Atom) = Byte4

  end do

end do
!----------------------------------------------------------------------*
! Print coordinates of the system                                      *
!----------------------------------------------------------------------*
write(6,*)
write(6,'(6X,A)') 'Cartesian coordinates in Angstrom:'
write(6,'(6X,A)') '-----------------------------------------------------'
write(6,'(6X,A)') 'No.  Label        X            Y            Z        '
write(6,'(6X,A)') '-----------------------------------------------------'
do iAt=1,iAll_Atom
  write(6,'(4X,I4,3X,A,2X,3F13.8)') iAt,AtomLbl(iAt),Angstr*W1(1:3,iAt)
end do
write(6,'(6X,A)') '-----------------------------------------------------'
write(6,'(6X,A,F14.8)') 'Nuclear repulsion energy =',PotNuc
call mma_deallocate(W1)

!----------------------------------------------------------------------*
! Normal exit                                                          *
!----------------------------------------------------------------------*
return

end subroutine PrCoor
