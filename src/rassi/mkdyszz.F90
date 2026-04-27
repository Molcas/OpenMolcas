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
! Copyright (C) 2018, Jesper Norell                                    *
!***********************************************************************

!****************************************************************
!  SUBROUTINE MKDYSZZ
!  PURPOSE: CALCULATE DYSON ORBITAL COEFFICIENTS FOR CI EXPANSIONS IN
!  BASIS FUNCTION BASE BASE Z,
!  IN ANALOGUE TO MKTDZZ FOR TRANSITION DENSITY MATRIX.
!****************************************************************
!  MODIFIED BY BRUNO TENORIO TO ADDRESS SYMMETRY
!  SEPTEMBER 2020
!****************************************************************

subroutine MKDYSZZ(CMOA,DYSAB,DYSZZ)

use Symmetry_Info, only: nIrrep
use rassi_data, only: NASH, NBASF, NCMO, NOSH
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: CMOA(NCMO), DYSAB(*)
real(kind=wp), intent(inout) :: DYSZZ(*)
integer(kind=iwp) :: BIOOFF, IBIO, IBIOFF, ISY1, IZZ, IZZOFF, NA1, NB1, NO1, SYMOFF
real(kind=wp) :: COEFF

! Re-express the DO coefficients in biorth basis DYSAB
! into atomic basis DYSZZ with help of CMOA that contains
! biorth orbitals in ZZ basis

SYMOFF = 0
IBIOFF = 0
IZZOFF = 0
do ISY1=1,nIrrep
  NO1 = NOSH(ISY1)
  NA1 = NASH(ISY1)
  NB1 = NBASF(ISY1)
  if (NA1 > 0) then
    do IBIO=1,NO1
      do IZZ=1,NB1
        BIOOFF = (IBIO-1)*NB1
        COEFF = DYSAB(IBIO+IBIOFF)*CMOA(SYMOFF+BIOOFF+IZZ)
        DYSZZ(IZZ+IZZOFF) = DYSZZ(IZZ+IZZOFF)+COEFF
      end do
    end do
  end if
  IZZOFF = NB1+IZZOFF
  IBIOFF = NO1+IBIOFF
  SYMOFF = (NO1*NB1)+SYMOFF
end do

end subroutine MKDYSZZ
