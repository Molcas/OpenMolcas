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

subroutine setsavvb_cvb(recn)

implicit real*8(a-h,o-z)
! ... Files/Hamiltonian available ...
logical, external :: tstfile_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
save recdef
data recdef/3200.2d0/

if (recn /= zero) return
do iadd=0,99
  if (.not. tstfile_cvb(recdef+dble(iadd))) then
    recn = recdef+dble(iadd)
    return
  end if
end do
recn = recdef+dble(99)

return

end subroutine setsavvb_cvb
