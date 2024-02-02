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

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: N, NM, L
real(kind=8), intent(in) :: W(N), T, dX, dY, dZ, H
real(kind=8), intent(out) :: MT, Z
complex(kind=8), intent(in) :: M(3,N,N)
integer :: i, j
real(kind=8) :: pB, dltw, S2, S1, mB, kB

! /// constants
mB = 0.4668643740_wp ! * in cm-1*T-1
kB = 0.69503560_wp   !   in cm^-1*K-1
pB = 0.0_wp
Z = 0.0_wp
MT = 0.0_wp
DLTW = 0.0_wp

do I=1,N
  pB = exp(-(W(I)-W(1))/kB/T)
  Z = Z+pB
  if (I <= NM) then
    ! case when I <= NM
    S2 = 0.0_wp
    S2 = dble(M(L,I,I))
    do J=NM+1,N
      DLTW = W(I)-W(J)
      S1 = 0.0_wp
      S1 = dble(M(L,I,J)*conjg(M(1,I,J)))*dX+dble(M(L,I,J)*conjg(M(2,I,J)))*dY+dble(M(L,I,J)*conjg(M(3,I,J)))*dZ
      S2 = S2-2.0_wp*mB*H*S1/DLTW
    end do ! J
  else ! I
    ! case when I > NM
    do J=1,N
      DLTW = W(I)-W(J)
      S1 = 0.0_wp
      S2 = 0.0_wp
      S1 = dble(M(L,I,J)*conjg(M(1,I,J)))*dX+dble(M(L,I,J)*conjg(M(2,I,J)))*dY+dble(M(L,I,J)*conjg(M(3,I,J)))*dZ

      if (abs(DLTW) < 1.D-3) then
        S2 = S2+1.0_wp*mB*H*S1/kB/T
      else
        S2 = S2-2.0_wp*mB*H*S1/DLTW
      end if
    end do ! J
  end if ! I
  MT = MT+pB*S2
end do !I

MT = MT/Z

return

end subroutine calcmagn2
