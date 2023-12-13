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

use Constants, only: Zero, One, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nparm, nnegeig
real(kind=wp), intent(in) :: heigval(nparm), gradp(nparm), hh, alfastart, alftol
real(kind=wp), intent(out) :: alfa
integer(kind=iwp) :: i
real(kind=wp) :: alfmax, alfmin, alfmx1, cnrm, cnrmax, cnrmin, olf, relfac

alfa = alfastart

! << Optimize alpha >>                             -1
! Norm of dX should be HH in:  dX = - (H - alpha I)   * G
olf = alfa
alfmin = alfa
alfmax = 1.0e2_wp+alfmin

relfac = One
do
  alfmin = olf
  alfmx1 = alfmax
  cnrmin = Zero
  cnrmax = Zero
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
  do
    alfa = Half*(alfmax+alfmin)
    cnrm = Zero
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
    if (abs(relfac*(alfmax-alfmin)) <= alftol) exit
  end do
  if (alfmax /= alfmx1) exit
  if (alfmax > 1.0e20_wp) then
    write(u6,*) ' Optimization of trust region size failed!'
    write(u6,*) ' Trust region size required :',hh
    write(u6,*) ' Min/max alpha values :',alfmin,alfmax
    write(u6,*) ' Min/max step sizes :',cnrmin,cnrmax
    call abend_cvb()
  end if
  alfmax = Ten*alfmax
  relfac = One/alfmax
end do
! Found the optimal alpha value, construct corresponding update:
alfa = Half*(alfmax+alfmin)

return

end subroutine optalf_cvb
