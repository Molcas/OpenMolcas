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

subroutine change3_cvb()

implicit real*8(a-h,o-z)
logical changed
! ... Change of dimensioning variables ...
logical, external :: chpcmp_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "rls_cvb.fh"
#include "WrkSpc.fh"

changed = .false.
! Spin functions coefficients (BIKCOF) + inverse (AIKCOF)
!  Get KBASISCVB if we don't know it already (eqv. to GETGUESS later):
!  Need to reserve enough space for both KBASIS and KBASISCVB
!  --> figure out which one needs most:
if (((kbasis > 2) .and. (kbasis /= 6)) .or. ((kbasiscvb > 2) .and. (kbasiscvb /= 6))) then
  kmost = 3
else if ((kbasis <= 2) .or. (kbasiscvb <= 2)) then
  kmost = 1
else
  kmost = 6
end if
if (chpcmp_cvb(kmost)) changed = .true.
if (changed) call touch_cvb('MEM3')

return

entry chop3_cvb()
if (release(3)) call mfreer_cvb(lb(1))
release(3) = .true.
release(4) = .false.

call icomb_cvb(nel,nbet,iretval1)
call icomb_cvb(nel,nbet-1,iretval2)
mxfns = iretval1-iretval2
if (kbasis == 5) call icomb_cvb(nel,nalf,mxfns)
call icomb_cvb(nel,nalf,mxdetvb)
if (((kbasis > 2) .and. (kbasis /= 6)) .or. ((kbasiscvb > 2) .and. (kbasiscvb /= 6))) then
  kmost = 3
else if ((kbasis <= 2) .or. (kbasiscvb <= 2)) then
  kmost = 1
else
  kmost = 6
end if
call bspset_cvb(kmost,1,need)
if (kmost == 3) then
  lb(1) = mstackr_cvb(1+need)
  lb(2) = mstackr_cvb(1+need)
else if (kmost == 1) then
  lb(1) = mstackr_cvb(1+need)
  lb(2) = lb(1)
else
  lb(1) = mstackr_cvb(1)
  lb(2) = lb(1)
end if
! Flag AIKCOF/BIKCOF as unset:
work(lb(1)) = zero
work(lb(2)) = zero

lb(3) = mstacki_cvb((nel+1)*(nel+1)*(nel+1))
lb(4) = mstacki_cvb((nel+1)*(nel+1))
lb(5) = mstacki_cvb((nel+1)*(nel+1))
lb(6) = mstacki_cvb((nel+1)*(nel+1))
call bspset_cvb(kbasiscvb,2,need)

return

end subroutine change3_cvb
