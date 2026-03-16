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

use Symmetry_Info, only: nSym => nIrrep
use Constants, only: Zero, One, Two
use rassi_data, only: NCMO, NBSQ, NBASF, NISH, NOSH

implicit none
real*8 CMO1(NCMO), CMO2(NCMO), DINAO(NBSQ)
integer ISTC, ISTD, ISYM, NI, NO, NB

DINAO(:) = Zero
ISTC = 1
ISTD = 1
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NO = NOSH(ISYM)
  NB = NBASF(ISYM)
  if (NI /= 0) call DGEMM_('N','T',NB,NB,NI,One,CMO1(ISTC),NB,CMO2(ISTC),NB,Zero,DINAO(ISTD),NB)
  ISTC = ISTC+NO*NB
  ISTD = ISTD+NB**2
end do
DINAO(:) = Two*DINAO(:)

end subroutine DIMAT
