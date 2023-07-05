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

subroutine ProjSym2(nAtoms,nCent,Ind,A,iDCRs,B,BqR,dB,dBqR)

use Slapaf_Info, only: jStab, nStab
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nAtoms, nCent, Ind(nCent), iDCRs(nCent)
real(kind=wp) :: A(3,nCent), B(3,nCent), BqR(3,nAtoms), dB(3,nCent,3,nCent), dBqR(3,nAtoms,3,nAtoms)
integer(kind=iwp) :: i, ixyz, j, jxyz
real(kind=wp) :: ATemp(3)
real(kind=wp), allocatable :: Tx(:,:)

#ifdef _DEBUGPRINT_
call RecPrt('B',' ',B,3,nCent)
call RecPrt('dB',' ',dB,3*nCent,3*nCent)
write(u6,*) iDCRs
#endif

! Set up the T-matrix

! Project away nonsymmetric displacements

call mma_allocate(Tx,3,nCent,Label='Tx')

call dcopy_(3*nCent,[One],0,Tx,1)
do i=1,nCent
  call NonSym(nStab(Ind(i)),jStab(0,Ind(i)),A(1,i),Tx(1,i))

  ! Rotate vector back to the unique center

  call OA(iDCRs(i),Tx(1:3,i),ATemp)
  Tx(:,i) = ATemp(:)
end do

! The T-matrix is now computed. Now create BqR and dBqR.

! Create BqR

call FZero(BqR,3*nAtoms)
do i=1,nCent
  do ixyz=1,3
    BqR(ixyz,Ind(i)) = BqR(ixyz,Ind(i))+Tx(ixyz,i)*B(ixyz,i)
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('BqR',' ',BqR,1,3*nAtoms)
#endif

! Create dBqR

call FZero(dBqR,(3*nAtoms)**2)
do i=1,nCent
  do ixyz=1,3

    do j=1,nCent
      do jxyz=1,3

        dBqR(ixyz,Ind(i),jxyz,Ind(j)) = dBqR(ixyz,Ind(i),jxyz,Ind(j))+Tx(ixyz,i)*dB(ixyz,i,jxyz,j)*Tx(jxyz,j)

      end do
    end do

  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('dBqR',' ',dBqR,3*nAtoms,3*nAtoms)
#endif

call mma_deallocate(Tx)

return

end subroutine ProjSym2
