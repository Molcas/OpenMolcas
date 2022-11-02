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
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Angstrom
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: iAll_atom, iAt, iAtom, iChAtom, iCo, iCoSet(0:7,0:7), iGen(3), iOper(0:7), iStab(0:7), MaxDCR, nAtoms, &
                     nCoSet, nGen, nStab, nSym
real(kind=wp) :: PotNuc
character(len=LenIn) :: Byte4
real(kind=wp), allocatable :: W1(:,:)
character(len=LenIn), allocatable :: AtomLbl(:)
integer(kind=iwp), external :: iChxyz

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
call mma_allocate(AtomLbl,MxAtom,label='AtomLbl')
call Get_cArray('Unique Atom Names',AtomLbl,LenIn*nAtoms)
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
  iChAtom = iChxyz(w1(:,iAtom),iGen,nGen)
  call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
  nCoSet = nSym/nStab
  Byte4 = AtomLbl(iAtom)

  do iCo=1,nCoSet-1

    iAll_Atom = iAll_Atom+1
    call OA(iCoSet(iCo,0),W1(:,iAtom),W1(:,iAll_Atom))
    AtomLbl(iAll_Atom) = Byte4

  end do

end do
!----------------------------------------------------------------------*
! Print coordinates of the system                                      *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,'(6X,A)') 'Cartesian coordinates in Angstrom:'
write(u6,'(6X,A)') '-----------------------------------------------------'
write(u6,'(6X,A)') 'No.  Label        X            Y            Z        '
write(u6,'(6X,A)') '-----------------------------------------------------'
do iAt=1,iAll_Atom
  write(u6,'(4X,I4,3X,A,2X,3F13.8)') iAt,AtomLbl(iAt),Angstrom*W1(:,iAt)
end do
write(u6,'(6X,A)') '-----------------------------------------------------'
write(u6,'(6X,A,F14.8)') 'Nuclear repulsion energy =',PotNuc
call mma_deallocate(AtomLbl)
call mma_deallocate(W1)

!----------------------------------------------------------------------*
! Normal exit                                                          *
!----------------------------------------------------------------------*
return

end subroutine PrCoor
