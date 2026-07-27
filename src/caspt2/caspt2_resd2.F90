!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CASPT2_ResD2(Mode,nRow,nCol,W1,W2,LDW,DIN,DIS)

use caspt2_global, only: imag_shift, real_shift, sigma_p_epsilon, sigma_p_exponent
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Mode, nRow, nCol, LDW
real(kind=wp), intent(inout) :: W1(LDW,nCol), W2(LDW,nCol)
real(kind=wp), intent(in) :: dIn(nRow), dIs(nCol)
integer(kind=iwp) :: i, j, p
real(kind=wp) :: delta, delta_inv, delta_ps, eps, expscal, scal, sigma

! See Eqs.(16), (17), and (22) in Stefano's sigma-regularization
! Eshift is the inverse of f(Delta;eps)
! real : Eshift = delta + epsilon
! imag : Eshift = delta + epsilon(^2)/delta
! sigma: Eshift = delta/(1-exp(-|delta/epsilon|^P))
! For sigma^P
! d(Eshift)/d(delta) = 1/(1-exp...) - delta/(1-...)
! 1st and 2nd terms are the ders of numerator and denominator
! The numerator   contribution is mode = 3
! The denominator contribution is mode = 2
! Note that Lagrangian is something like (scaling factors can be
! case-dependent)
! 2<1|V|0> + <1|H0-E0|1> + <lambda|V|0> + <lambda|H0-E0+Eshift|1>
! The correlated density comes from the 2nd and 4th terms
eps = Zero
if (sigma_p_epsilon /= Zero) eps = sigma_p_epsilon

do j=1,nCol
  do i=1,nRow
    scal = Zero
    select case (mode)
      case (1)
        ! energy denominator plus real shift
        delta = dIn(i)+dIs(j)+real_shift
        ! inverse denominator plus imaginary shift
        delta_inv = delta/(delta**2+imag_shift**2)
        ! multiply by (inverse) sigma-p regularizer
        eps = sigma_p_epsilon
        p = sigma_p_exponent
        if (eps > Zero) then
          sigma = One/eps**p
          delta_inv = delta_inv*(One-exp(-sigma*abs(delta)**p))
        end if
            !! The following SCAL is the actual residual
        scal = One-(dIn(i)+dIs(j))*delta_inv
            !! Another scaling is required for lambda
        scal = -scal*delta_inv

        W1(i,j) = scal*W1(i,j)
        W2(i,j) = scal*W2(i,j)
      case (2)
        if (imag_shift /= Zero) then
          scal = imag_shift/(dIn(i)+dIs(j))
          W1(i,j) = scal*W1(i,j)
          W2(i,j) = scal*W2(i,j)
        else if (eps /= Zero) then
          ! derivative of the denominator of sigma-p CASPT2
          ! always real_shift = imag_shift = 0
          delta = dIn(i)+dIs(j)
          ! multiply by (inverse) sigma-p regularizer
          p = sigma_p_exponent
          sigma = One/eps**p
          delta_ps = (delta**p)*sigma
          expscal = exp(-abs(delta_ps))
          delta_inv = One/(One-expscal)
          ! Correct only for diagonal elements
          ! Use the canonical constraint for off-diagonal elements
          W1(i,j) = delta_inv*W1(i,j)*expscal*delta/eps*abs(delta/eps)**real(p-1,kind=wp)*real(p,kind=wp)
          W2(i,j) = delta_inv*W2(i,j)
        end if
      case (3)
        ! derivative of the numerator of sigma-p CASPT2
        ! called only for sigma-p CASPT2
        ! always real_shift = imag_shift = 0
        delta = dIn(i)+dIs(j)
        ! multiply by (inverse) sigma-p regularizer
        p = sigma_p_exponent
        sigma = One/eps**p
        expscal = exp(-sigma*abs(delta)**p)
        delta_inv = One/(One-expscal)
        ! The relevant terms are:
        ! <1|H0-E0|1> + <lambda|H0-E0+Eshift|1>
        ! This subroutine modifies the lambda for <lambda|Eshift|1>
        ! so scaling only the lambda is correct.
        W1(i,j) = delta_inv*W1(i,j)
    end select
  end do
end do

end subroutine CASPT2_ResD2
