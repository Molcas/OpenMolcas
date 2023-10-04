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

!IFG trivial nel, kbasiscvb, ip
subroutine mkbiks_cvb()

use casvb_global, only: aikcof, bikcof

implicit none
#include "main_cvb.fh"
#include "print_cvb.fh"

call biks_cvb(aikcof,bikcof,nel,kbasiscvb,associated(bikcof,aikcof),ip(1))

return

end subroutine mkbiks_cvb
