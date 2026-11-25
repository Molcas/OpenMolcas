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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine H0DIAG_CASPT2(ISYCI,DIAG,nDiag)
! PURPOSE: FORM AN ARRAY OF DIAGONAL HAMILTONIAN MATRIX ELEMENTS
! FOR THE SPECIFIED TOTAL SYMMETRY ISYCI: sum_lev n_lev ETA(lev) over the
! CSFs (CI interface, CI_Occ_Diag)

use CI_Lib_States, only: CI_Occ_Diag
use caspt2_module, only: ETA
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ISYCI, nDiag
real(kind=wp), intent(out) :: DIAG(nDiag)
integer(kind=iwp), parameter :: istate = 1

call CI_Occ_Diag(istate,ISYCI,ETA,DIAG,nDiag)

end subroutine H0DIAG_CASPT2
