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

logical function ifcasci_cvb()

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "applyh_cvb.fh"

call f_inquire('JOBOLD',ifcasci_cvb)
! In variational calculations, CI vectors will be of no use
! unless they correspond to MOs identical to the present ones:
if (variat .and. ((invec_cvb /= 3) .or. (nmcscf > 1))) ifcasci_cvb = .false.

return

end function ifcasci_cvb
