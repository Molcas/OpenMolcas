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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CLagDX_Sigma(iSym,iCase,nIN,nIS,nAS,jSta,nCol,VEC1,VEC2,EIG,W1)

use EQSOLV, only: IDBMAT
use caspt2_global, only: LUSBT, sigma_p_epsilon, sigma_p_exponent
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym, iCase, nIN, nIS, nAS, jSta, nCol
real(kind=wp), intent(in) :: VEC1(nIN,nCol), VEC2(nIN,nCol), EIG(nIN)
real(kind=wp), intent(inout) :: W1(nIN,nIN)
integer(kind=iwp) :: idB, iICB, iIS, jICB
real(kind=wp) :: delta, Dummy(1), expscal, fact
real(kind=wp), allocatable :: derscal(:), DIS(:), scal(:)
real(kind=wp), parameter :: ThrDeg = 1.0e-8_wp

! Derivative of the sigma^p shift with respect to the B matrix

call mma_allocate(DIS,nIS,Label='DIS')
call mma_allocate(scal,nIN,Label='scal')
call mma_allocate(derscal,nIN,Label='derscal')

idB = IDBMAT(iSym,iCase)
call DDAFILE(LUSBT,0,Dummy,nAS,idB)
call DDAFILE(LUSBT,2,DIS,nIS,idB)

do iIS=1,nCol
  do iICB=1,nIN
    delta = EIG(iICB)+DIS(jSta+iIS-1)
    expscal = exp(-abs(delta/sigma_p_epsilon)**sigma_p_exponent)
    scal(iICB) = delta*expscal/(One-expscal)
    derscal(iICB) = (expscal-sigma_p_exponent*(abs(delta/sigma_p_epsilon)**sigma_p_exponent)*expscal/(One-expscal))/(one-expscal)
  end do
  do jICB=1,nIN
    do iICB=1,nIN
      if (abs(EIG(iICB)-EIG(jICB)) > ThrDeg) then
        fact = (scal(iICB)-scal(jICB))/(EIG(iICB)-EIG(jICB))
      else
        fact = Half*(derscal(iICB)+derscal(jICB))
      end if
      W1(iICB,jICB) = W1(iICB,jICB)+Half*(VEC2(iICB,iIS)*VEC1(jICB,iIS)+VEC1(iICB,iIS)*VEC2(jICB,iIS))*fact
    end do
  end do
end do

call mma_deallocate(DIS)
call mma_deallocate(scal)
call mma_deallocate(derscal)

end subroutine CLagDX_Sigma
