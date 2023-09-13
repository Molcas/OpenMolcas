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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine optalf_cvb(heigval,gradp,nparm,hh,alfa,nnegeig,alfastart,alftol)

implicit real*8(a-h,o-z)
dimension heigval(nparm), gradp(nparm)
save zero, half, one
data zero/0.d0/,half/0.5d0/,one/1.d0/

alfa = alfastart

! << Optimize alpha >>                             -1
! Norm of dX should be HH in:  dX = - (H - alpha I)   * G
olf = alfa
alfmin = alfa
alfmax = 1.d2+alfmin

relfac = one
700 alfmin = olf
alfmx1 = alfmax
cnrmin = zero
cnrmax = zero
do i=1,nnegeig
  cnrmin = cnrmin+(gradp(i)/(heigval(i)-alfmin))**2
  cnrmax = cnrmax+(gradp(i)/(heigval(i)-alfmax))**2
end do
do i=nnegeig+1,nparm
  cnrmin = cnrmin+(gradp(i)/(heigval(i)+alfmin))**2
  cnrmax = cnrmax+(gradp(i)/(heigval(i)+alfmax))**2
end do
cnrmin = sqrt(cnrmin)
cnrmax = sqrt(cnrmax)
900 alfa = half*(alfmax+alfmin)
cnrm = zero
do i=1,nnegeig
  cnrm = cnrm+(gradp(i)/(heigval(i)-alfa))**2
end do
do i=nnegeig+1,nparm
  cnrm = cnrm+(gradp(i)/(heigval(i)+alfa))**2
end do
cnrm = sqrt(cnrm)
if (cnrm < hh) then
  alfmax = alfa
  cnrmax = cnrm
else
  alfmin = alfa
  cnrmin = cnrm
end if
if (abs(relfac*(alfmax-alfmin)) > alftol) goto 900
if (alfmax == alfmx1) then
  if (alfmax > 1d20) then
    write(6,*) ' Optimization of trust region size failed!'
    write(6,*) ' Trust region size required :',hh
    write(6,*) ' Min/max alpha values :',alfmin,alfmax
    write(6,*) ' Min/max step sizes :',cnrmin,cnrmax
    call abend_cvb()
  end if
  alfmax = 1.d1*alfmax
  relfac = one/alfmax
  goto 700
end if
! Found the optimal alpha value, construct corresponding update:
alfa = half*(alfmax+alfmin)

return

end subroutine optalf_cvb
