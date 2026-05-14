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

use Constants, only: Zero, Two, cm_s, hPlanck, kBoltzmann, mBohr
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, NM, L
real(kind=wp), intent(in) :: W(N), T, H, dX, dY, dZ
complex(kind=wp), intent(in) :: M(3,N,N)
real(kind=wp), intent(out) :: MT, Z
integer(kind=iwp) :: i, j
real(kind=wp) :: dltw, pB, S1, S2
real(kind=wp), parameter :: kB = kBoltzmann/(cm_s*hPlanck), & ! in cm-1*K-1
                            mB = mBohr/(cm_s*hPlanck) ! in cm-1*T-1

Z = Zero
MT = Zero

do I=1,N
  pB = exp(-(W(I)-W(1))/kB/T)
  Z = Z+pB
  if (I <= NM) then
    ! case when I <= NM
    S2 = real(M(L,I,I))
    do J=NM+1,N
      DLTW = W(I)-W(J)
      S1 = real(M(L,I,J)*conjg(M(1,I,J)))*dX+real(M(L,I,J)*conjg(M(2,I,J)))*dY+real(M(L,I,J)*conjg(M(3,I,J)))*dZ
      S2 = S2-Two*mB*H*S1/DLTW
    end do ! J
  else ! I
    ! case when I > NM
    do J=1,N
      DLTW = W(I)-W(J)
      S1 = real(M(L,I,J)*conjg(M(1,I,J)))*dX+real(M(L,I,J)*conjg(M(2,I,J)))*dY+real(M(L,I,J)*conjg(M(3,I,J)))*dZ
      S2 = Zero

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
