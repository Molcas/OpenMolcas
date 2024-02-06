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

subroutine calcmagn2(N,NM,W,T,H,M,dX,dY,dZ,L,MT,Z)
! this subroutine sums the contribution to molar magnetization from
! different states:
!   -- from the states obtained in exact Zeeman diagonalization (NM) and
!   -- from the states higher in energy via second order perturbation
!----------------------------------------------------------------------
! Input data:
!    N  -  total number of states
!   NM  -  number of states entering the exact Zeeman diagonalization
!         (NM <= N)
!  W(N) -  energy of the input states

use Constants, only: Zero, Two, cLight, kBoltzmann, mBohr, rPlanck
use Definitions, only: wp

implicit none
integer, intent(in) :: N, NM, L
real(kind=8), intent(in) :: W(N), T, dX, dY, dZ, H
real(kind=8), intent(out) :: MT, Z
complex(kind=8), intent(in) :: M(3,N,N)
integer :: i, j
real(kind=8) :: pB, dltw, S2, S1
real(kind=wp), parameter :: kB = kBoltzmann/(cLight*rPlanck*1.0e2_wp), & ! in cm-1*K-1
                            mB = mBohr/(cLight*rPlanck*1.0e2_wp) ! in cm-1*T-1

pB = Zero
Z = Zero
MT = Zero
DLTW = Zero

do I=1,N
  pB = exp(-(W(I)-W(1))/kB/T)
  Z = Z+pB
  if (I <= NM) then
    ! case when I <= NM
    S2 = Zero
    S2 = real(M(L,I,I))
    do J=NM+1,N
      DLTW = W(I)-W(J)
      S1 = Zero
      S1 = real(M(L,I,J)*conjg(M(1,I,J)))*dX+real(M(L,I,J)*conjg(M(2,I,J)))*dY+real(M(L,I,J)*conjg(M(3,I,J)))*dZ
      S2 = S2-Two*mB*H*S1/DLTW
    end do ! J
  else ! I
    ! case when I > NM
    do J=1,N
      DLTW = W(I)-W(J)
      S1 = Zero
      S2 = Zero
      S1 = real(M(L,I,J)*conjg(M(1,I,J)))*dX+real(M(L,I,J)*conjg(M(2,I,J)))*dY+real(M(L,I,J)*conjg(M(3,I,J)))*dZ

      if (abs(DLTW) < 1.0e-3_wp) then
        S2 = S2+mB*H*S1/kB/T
      else
        S2 = S2-Two*mB*H*S1/DLTW
      end if
    end do ! J
  end if ! I
  MT = MT+pB*S2
end do ! I

MT = MT/Z

return

end subroutine calcmagn2
