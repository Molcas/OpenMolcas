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
!*************************************************************
!** Routines for imposing constraints on VB wfn. parameters **
!*************************************************************
!*********************
!** Set-up routines **
!*********************

subroutine mkconstruc_cvb()

use casvb_global, only: iconstruc, ipermzeta, izeta, orbs, symelm, tconstr

implicit none

call setipermzeta_cvb(ipermzeta,orbs,symelm,izeta)
if (iconstruc == 2) call construc2_cvb(tconstr)

return

end subroutine mkconstruc_cvb
