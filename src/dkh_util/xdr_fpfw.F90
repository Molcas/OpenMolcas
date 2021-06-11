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

subroutine XDR_fpFW(n,s,t,v,w,Tr,Bk,EL,ES,OL,OS,P,E0,A,B,R,clight)
! Transform to moment space ( eigenfunction space of non-relativistic kinetic matrix T )
!   and apply the free-particle Foldy-Wothuysen transformation to the Fock matrices

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Two
use Definitions, only: wp, iwp

implicit none
! Input
!   w   aka pVp
! Output :
!   Tr     transform matrix to moment space
!   Bk     back transform matrix of moment space =Tr^{-1}
!   ( EL OL )
!   ( OS ES ) represent the structure of potential matrix in fpFW space
!   P      kinetic matrix (diagonal) in fpFW space, aka Ep
!   E0     =P-mc^{2}
!   A,B,R  kinetic factors (used to construct the fpFW transformation in moment space)
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: s(n,n), t(n,n), v(n,n), w(n,n), clight
real(kind=wp), intent(out) :: Tr(n,n), Bk(n,n), EL(n,n), ES(n,n), OL(n,n), OS(n,n), P(n), E0(n), A(n), B(n), R(n)
integer(kind=iwp) :: i, j, lwork, info
real(kind=wp) :: av, aw
real(kind=wp), allocatable :: tmp(:), sW(:), sA(:,:), sB(:,:)

! Diagonalization

lwork = 8*n
call mma_allocate(tmp,lwork,label='Tmp')
call mma_allocate(sW,n,label='Eig')
Tr(:,:) = t(:,:)
Bk(:,:) = s(:,:)
call dsygv_(1,'V','L',n,Tr,n,Bk,n,sW,tmp,lwork,info)

! Transform potential matrices to moment space

call mma_allocate(sA,n,n,label='TmpA')
call mma_allocate(sB,n,n,label='TmpB')
call dmxma(n,'C','N',Tr,v,Bk,One)
call dmxma(n,'N','N',Bk,Tr,sA,One)
call dmxma(n,'C','N',Tr,w,Bk,One)
call dmxma(n,'N','N',Bk,Tr,sB,One)

! Calculate kinetic moment factors

do i=1,n
  P(i) = sqrt(clight*clight+sW(i)*Two)*clight
  A(i) = sqrt((P(i)+clight*clight)/(P(i)+P(i)))
  B(i) = clight/sqrt(Two*P(i)*(P(i)+clight*clight))
  R(i) = sqrt(sW(i)*Two)*clight/(P(i)+clight*clight)
  ! E0 equal to Ep-c^{2}, but this formula is more stable, especially E0 close to zero
  E0(i) = Two*sW(i)*clight*clight/(P(i)+clight*clight)
end do

! Multiply potential matrices with moment factors
!   aka the fpFW transformation in moment space

do i=1,n
  do j=1,n
    av = sA(j,i)*A(i)*A(j)
    aw = sB(j,i)*B(i)*B(j)
    EL(j,i) = av+aw
    ES(j,i) = aw/R(i)/R(j)+av*R(i)*R(j)
    OL(j,i) = aw/R(i)-av*R(i)
    OS(j,i) = aw/R(j)-av*R(j)
  end do
end do
Bk(:,:) = Tr(:,:)

! Make back transform matrix

call XDR_dmatinv(Bk,n)

! Free temp memories

call mma_deallocate(tmp)
call mma_deallocate(sW)
call mma_deallocate(sA)
call mma_deallocate(sB)

return

end subroutine XDR_fpFW
