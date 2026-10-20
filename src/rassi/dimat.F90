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
! Copyright (C) 2000, Per Ake Malmqvist                                *
!***********************************************************************
!****************************************************************
!  PROGRAM RASSI        PER-AAKE MALMQVIST 2000-06-30
!  SUBROUTINE DIMAT
!  CONSTRUCT A DENSITY MATRIX DINAO FOR THE INACTIVE ORBITALS.
!  IT IS RETURNED IN SYMMETRY-BLOCKED SQUARED FORMAT.
!****************************************************************

subroutine DIMAT(CMO1,CMO2,DINAO)

use rassi_data, only: NCMO, NBSQ, NBASF, NISH, NOSH
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: CMO1(NCMO), CMO2(NCMO)
real(kind=wp), intent(out) :: DINAO(NBSQ)
integer(kind=iwp) :: ISTC, ISTD, ISYM, NB, NI, NO

ISTC = 1
ISTD = 1
do ISYM=1,nIrrep
  NI = NISH(ISYM)
  NO = NOSH(ISYM)
  NB = NBASF(ISYM)
  if (NI == 0) then
    DINAO(ISTD:ISTD+NB**2-1) = Zero
  else
    call DGEMM_('N','T',NB,NB,NI,Two,CMO1(ISTC),NB,CMO2(ISTC),NB,Zero,DINAO(ISTD),NB)
  end if
  ISTC = ISTC+NO*NB
  ISTD = ISTD+NB**2
end do

end subroutine DIMAT
