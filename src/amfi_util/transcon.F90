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

subroutine transcon(contold,idim1,idim2,ovlp,contnew,nprim,ncont)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idim1, idim2, nprim, ncont
real(kind=wp), intent(in) :: contold(idim1,idim2), ovlp(idim1,idim1)
real(kind=wp), intent(out) :: contnew(nprim,ncont)
integer(kind=iwp) :: ICONT, Irun, Jrun
real(kind=wp) :: xnorm

!write(u6,*) 'begin transcon nprim,ncont ',nprim,ncont
!bs copy old contraction coefficients in dense form
contnew(:,:) = contold(1:nprim,1:ncont)
!bs ensure normalization
do ICONT=1,ncont
  xnorm = Zero
  do Jrun=1,nprim
    do Irun=1,nprim
      xnorm = xnorm+contnew(Irun,ICONT)*contnew(Jrun,ICONT)*ovlp(Irun,Jrun)
      !write(u6,*) 'Icont,jrun,irun,xnorm ',icont,jrun,irun,xnorm
    end do
  end do
  !write(u6,*) 'ICONT ',ICONT,xnorm
  xnorm = One/sqrt(xnorm)
  !bs scale with normalization factor
  contnew(:,ICONT) = xnorm*contnew(:,ICONT)
end do
!write(u6,*) 'end transcon nprim,ncont ',nprim,ncont

return

end subroutine transcon
