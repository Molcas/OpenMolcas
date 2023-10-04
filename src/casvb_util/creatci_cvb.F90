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

subroutine creatci_cvb(inumber)

use casvb_global, only: civbvecs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: inumber
#include "main_cvb.fh"
logical(kind=iwp), parameter :: debug = .false.

civbvecs(1,inumber) = real(inumber,kind=wp)
iform_ci(inumber) = nint(civbvecs(0,inumber))
if (debug) then
  write(u6,*) ' Creating CI vector :',inumber
  write(u6,*) ' Format             :',iform_ci(inumber)
end if

return

end subroutine creatci_cvb
