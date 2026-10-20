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

subroutine zorder(ndimen,ldv,vecre,vecim,arrre,switch)
! Order the eigenvalues in increasing sequence:
! Note that special care is taken to define some order in case
! of degeneracy.
! switch controls whether we assume an array or vector
! of eigenvalues

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIMEN, LDV, switch
real(kind=wp), intent(inout) :: VECRE(LDV,*), VECIM(LDV,*), ARRRE(NDIMEN,*)
integer(kind=iwp) :: I, ISEL, J, K
real(kind=wp) :: EDIFF, ESEL, EVAL, O_i, O_j, VIKI, VIKJ, VRKI, VRKJ
real(kind=wp), parameter :: Thr_EDiff = 1.0e-10_wp

do I=1,NDIMEN-1
  ESEL = ARRRE(I,I**(switch))
  ISEL = I

  O_i = Zero
  do K=1,LDV
    O_i = O_i+real(K,kind=wp)*(VECRE(K,I)**2+VECIM(K,I)**2)
  end do

  ! Check against other eigenvalues

  do J=I+1,NDIMEN
    EVAL = ARRRE(J,J**switch)

    EDIFF = abs(EVAL-ESEL)
    if ((EVAL < ESEL) .and. (EDIFF > Thr_EDiff)) then
      ISEL = J
      ESEL = EVAL
    else if (EDIFF < Thr_EDiff) then
      O_j = Zero
      do K=1,LDV
        O_j = O_j+real(K,kind=wp)*(VECRE(K,J)**2+VECIM(K,J)**2)
      end do
      if (O_j > O_i) then
        ISEL = J
        ESEL = EVAL
      end if
    end if
  end do

  ! Swap here!

  if (ISEL /= I) then
    do K=1,LDV
      VRKI = VECRE(K,I)
      VRKJ = VECRE(K,ISEL)
      VIKI = VECIM(K,I)
      VIKJ = VECIM(K,ISEL)
      VECRE(K,I) = VRKJ
      VECRE(K,ISEL) = VRKI
      VECIM(K,I) = VIKJ
      VECIM(K,ISEL) = VIKI
    end do
    ARRRE(ISEL,ISEL**switch) = ARRRE(I,I**switch)
    ARRRE(I,I**switch) = ESEL
  end if
end do

end subroutine zorder
