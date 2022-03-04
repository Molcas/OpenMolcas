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

implicit real*8(a-h,o-z)
dimension contold(idim1,idim2), contnew(nprim,ncont),ovlp(idim1,idim1)

!write(6,*) 'begin transcon nprim,ncont ',nprim,ncont
!bs copy old contraction coefficients in dense form to common block
do Jrun=1,ncont
  do Irun=1,nprim
    contnew(Irun,Jrun) = contold(Irun,Jrun)
  end do
end do
!bs ensure normalization
do ICONT=1,ncont
  xnorm = 0d0
  do Jrun=1,nprim
    do Irun=1,nprim
      xnorm = xnorm+contnew(Irun,ICONT)*contnew(Jrun,ICONT)*ovlp(Irun,Jrun)
      !write(6,*) 'Icont,jrun,irun,xnorm ',icont,jrun,irun,xnorm
    end do
  end do
  !write(6,*) 'ICONT ',ICONT,xnorm
  xnorm = 1d0/sqrt(xnorm)
  !bs scale with normalization factor
  do Irun=1,nprim
    contnew(Irun,ICONT) = xnorm*contnew(Irun,ICONT)
  end do
end do
!write(6,*) 'end transcon nprim,ncont ',nprim,ncont

return

end subroutine transcon
