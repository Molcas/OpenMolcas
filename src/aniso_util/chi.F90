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

subroutine chi(M1,M2,E,N,T,Z,X)
! computes the chi  of one site
! definition of the variables:
!     N -- number of states included in the calculation of chi, Integer, input
!    M1 -- moment 1, Complex*16, (3,N,N) array, input
!    M2 -- moment 2, Complex*16, (3,N,N) array, input
!     E -- energy of the N states, Real(kind=wp) ::, (N) array, input
!     T -- temperature at which the chi is computed, Real(kind=wp) ::, input
!     Z -- statistical sum according to Boltzmann distribution law, Real(kind=wp) ::, output
!     X -- susceptibility tensor, Real(kind=wp) ::, (3,3) array, output
!--------
! temporary (local) variables:
! iS,jS -- denote states over which the chi is computed
!    pB -- partial Boltzmann population of a given state, Real(kind=wp) ::
!    dE -- energy diference E(i)-E(j)

use Constants, only: Zero, One, Two, cm_s, hPlanck, kBoltzmann
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
complex(kind=wp), intent(in) :: M1(3,N,N), M2(3,N,N)
real(kind=wp), intent(in) :: E(N), T
real(kind=wp), intent(out) :: Z, X(3,3)
integer(kind=iwp) :: i, iS, j, jS
real(kind=wp) :: c2(3,3), dE, F, pB, R
real(kind=wp), parameter :: boltz_k = kBoltzmann/(cm_s*hPlanck) ! in cm-1*K-1

Z = Zero
X(:,:) = Zero

do iS=1,N
  ! first loop over all states
  c2(:,:) = Zero
  ! pB = statistical sum for state iS at temperature T
  pB = exp(-E(iS)/boltz_k/T)
  ! accumulate the total statistical sum Z
  Z = Z+pB
  do jS=1,N
    ! second loop over all states
    dE = E(iS)-E(jS)
    ! set the multiplication factor:
    if (abs(dE) < 1.0e-3_wp) then
      F = One
    else
      F = -Two*boltz_k*T/dE
    end if
    ! accumulate the contributions to the X tensor:
    do i=1,3
      do j=1,3
        R = real(M1(i,iS,jS)*conjg(M2(j,iS,jS)))
        c2(i,j) = c2(i,j)+F*R
      end do ! j
    end do ! i
  end do ! jS

  ! add the (iS) contribution to the X tensor, according to its
  ! Boltzmann distribution:
  X(:,:) = X(:,:)+pB*c2(:,:)
end do ! iS

! scale the total tensor by the total statistical sum Z:
X(:,:) = X(:,:)/Z

return

end subroutine chi
