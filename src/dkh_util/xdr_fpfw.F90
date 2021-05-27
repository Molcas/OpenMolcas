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
integer(kind=iwp) :: i, j, k, itmp, lwork, iW, iA, iB, info
real(kind=wp) :: av, aw
#include "WrkSpc.fh"

! Diagonalization

lwork = 8*n
call getmem('Tmp  ','ALLOC','REAL',itmp,lwork)
call getmem('Eig  ','ALLOC','REAL',iW,n+4)
k = 0
do i=1,n
  do j=1,n
    Tr(j,i) = t(j,i)
    Bk(j,i) = s(j,i)
    k = k+1
  end do
end do
call dsygv_(1,'V','L',n,Tr,n,Bk,n,Work(iW),Work(itmp),lwork,info)

! Transform potential matrices to moment space

call getmem('TmpA ','ALLOC','REAL',iA,n*n+4)
call getmem('TmpB ','ALLOC','REAL',iB,n*n+4)
call dmxma(n,'C','N',Tr,v,Bk,One)
call dmxma(n,'N','N',Bk,Tr,Work(iA),One)
call dmxma(n,'C','N',Tr,w,Bk,One)
call dmxma(n,'N','N',Bk,Tr,Work(iB),One)

! Calculate kinetic moment factors

do i=1,n
  P(i) = sqrt(clight*clight+Work(iW+i-1)*Two)*clight
  A(i) = sqrt((P(i)+clight*clight)/(P(i)+P(i)))
  B(i) = clight/sqrt(Two*P(i)*(P(i)+clight*clight))
  R(i) = sqrt(Work(iW+i-1)*Two)*clight/(P(i)+clight*clight)
  ! E0 equal to Ep-c^{2}, but this formula is more stable, especially E0 close to zero
  E0(i) = Two*Work(iW+i-1)*clight*clight/(P(i)+clight*clight)
end do

! Multiply potential matrices with moment factors
!   aka the fpFW transformation in moment space

k = 0
do i=1,n
  do j=1,n
    av = Work(iA+k)*A(i)*A(j)
    aw = Work(iB+k)*B(i)*B(j)
    EL(j,i) = av+aw
    ES(j,i) = aw/R(i)/R(j)+av*R(i)*R(j)
    OL(j,i) = aw/R(i)-av*R(i)
    OS(j,i) = aw/R(j)-av*R(j)
    Bk(j,i) = Tr(j,i)
    k = k+1
  end do
end do

! Make back transform matrix

call XDR_dmatinv(Bk,n)

! Free temp memories

call getmem('Tmp  ','FREE','REAL',itmp,lwork)
call getmem('Eig  ','FREE','REAL',iW,n+4)
call getmem('TmpA ','FREE','REAL',iA,n*n+4)
call getmem('TmpB ','FREE','REAL',iB,n*n+4)

return

end subroutine XDR_fpFW
