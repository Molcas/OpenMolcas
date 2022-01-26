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

subroutine DIAG_DKH(A,N,EIG,EW,SINV,AUX,IC)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: N, IC
real(kind=wp) :: A(N*(N+1)/2), EIG(N,N), EW(N), SINV(N,N), AUX(N,N)
integer(kind=iwp) :: I, IJ, J, K, L
real(kind=wp) :: TOL
#ifdef MOLPRO
!MR ATTENTION, THE SCRATCH ARRAY TMP IS NOT PROPERLY ALLOCATED
real(kind=wp) :: TMP(6*N)
#endif

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUX(I,J) = A(IJ)
    AUX(J,I) = A(IJ)
  end do
end do
do K=1,N
  do J=1,N
    EIG(K,J) = Zero
    do L=1,J
      EIG(K,J) = EIG(K,J)+AUX(K,L)*SINV(L,J)
    end do
  end do
end do
do I=1,N
  do J=1,I
    AUX(I,J) = Zero
    do K=1,I
      AUX(I,J) = AUX(I,J)+SINV(K,I)*EIG(K,J)
    end do
    AUX(J,I) = AUX(I,J)
  end do
end do

TOL = 1.0e-80_wp
#ifdef MOLPRO
do I=1,N
  do J=1,N
    EIG(J,I) = AUX(J,I)
  end do
end do
call diag2(n,n,ew,eig)
!call dsyev_('V','L',N,EIG,N,EW,TMP,6*N,INFO)
#else
call JACOB_REL(AUX,EIG,EW,N,TOL,IC)
#endif

return

end subroutine DIAG_DKH
