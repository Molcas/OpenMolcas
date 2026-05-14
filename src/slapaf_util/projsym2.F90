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
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nCent, Ind(nCent), iDCRs(nCent)
real(kind=wp), intent(in) :: A(3,nCent), B(3,nCent), dB(3,nCent,3,nCent)
real(kind=wp), intent(out) :: BqR(3,nAtoms), dBqR(3,nAtoms,3,nAtoms)
integer(kind=iwp) :: i, j, jxyz
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

Tx(:,:) = One
do i=1,nCent
  call NonSym(nStab(Ind(i)),jStab(0,Ind(i)),A(:,i),Tx(:,i))

  ! Rotate vector back to the unique center

  call OA(iDCRs(i),Tx(1:3,i),ATemp)
  Tx(:,i) = ATemp(:)
end do

! The T-matrix is now computed. Now create BqR and dBqR.

! Create BqR

BqR(:,:) = Zero
do i=1,nCent
  BqR(:,Ind(i)) = BqR(:,Ind(i))+Tx(:,i)*B(:,i)
end do
#ifdef _DEBUGPRINT_
call RecPrt('BqR',' ',BqR,1,3*nAtoms)
#endif

! Create dBqR

dBqR(:,:,:,:) = Zero
do j=1,nCent
  do jxyz=1,3

    do i=1,nCent

      dBqR(:,Ind(i),jxyz,Ind(j)) = dBqR(:,Ind(i),jxyz,Ind(j))+Tx(:,i)*dB(:,i,jxyz,j)*Tx(jxyz,j)

    end do

  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('dBqR',' ',dBqR,3*nAtoms,3*nAtoms)
#endif

call mma_deallocate(Tx)

return

end subroutine ProjSym2
