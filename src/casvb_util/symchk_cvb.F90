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

subroutine symchk_cvb()

use Definitions, only: iwp

implicit none
logical(kind=iwp), external :: recinpcmp_cvb, up2date_cvb ! ... Make: up to date? ...

if (up2date_cvb('SYMINIT')) then
  ! iorts:
  if (recinpcmp_cvb(14)) call touch_cvb('ORBFREE')
  ! irots:
  if (recinpcmp_cvb(15)) call touch_cvb('ORBFREE')
  ! iorbrel:
  if (recinpcmp_cvb(9)) then
    call touch_cvb('SYMINIT')
    call touch_cvb('ORBFREE')
  end if
  ! ifxorb:
  if (recinpcmp_cvb(11)) then
    call touch_cvb('SYMINIT')
    call touch_cvb('ORBFREE')
  end if
end if
if (up2date_cvb('CONSTRUC')) then
  ! ifxstr:
  if (recinpcmp_cvb(12)) then
    call touch_cvb('CONSTRUC')
    call touch_cvb('CIFREE')
  end if
  ! idelstr:
  if (recinpcmp_cvb(13)) then
    call touch_cvb('CONSTRUC')
    call touch_cvb('CIFREE')
  end if
  ! izeta:
  if (recinpcmp_cvb(16)) then
    call touch_cvb('CONSTRUC')
    call touch_cvb('CIFREE')
  end if
end if

return

end subroutine symchk_cvb
