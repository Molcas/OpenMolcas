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
! Copyright (C) 1989, Per Ake Malmqvist                                *
!***********************************************************************
!****************************************************************
!  PROGRAM RASSI        PER-AAKE MALMQVIST
!  SUBROUTINE FINDT     IBM-3090 RELEASE 89 01 31
!  GIVEN CMO1 AND CMO2, FIND SUITABLE COUEFFICIENTS FOR THE
!  SINGLE-ORBITAL TRANSFORMATION SEQUENCE FOR EACH SET, RETURNED
!  IN ARRAYS TRA1 AND TRA2, FOR TRANSFORMATION INTO BIORTHONORMAL
!  ORBITALS.
!****************************************************************

subroutine FINDT(CMO1,CMO2,TRA1,TRA2)

use Symmetry_Info, only: nIrrep
use rassi_data, only: NBASF, NCMO, NOSH, NSXY, NTRA
use Cntrl, only: PRORB, PRSXY, PRTRA
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMO1(NCMO), CMO2(NCMO), TRA1(NTRA), TRA2(NTRA)
integer(kind=iwp) :: NCXA, NCYB
real(kind=wp), allocatable :: CXA(:), CYB(:), SXY(:)

! BEGIN BY COMPUTING MO OVERLAPS, SXY.
! MO OVERLAP MATRIX HAS SAME STRUCTURE AS TRA1,TRA2:
NSXY = NTRA
call mma_allocate(SXY,NSXY,Label='SXY')
call MKSXY(CMO1,CMO2,SXY)
if (PRSXY) then
  write(u6,*)
  call WRMAT('MO OVERLAP MATRIX IN ORIGINAL ORBITAL BASES:',1,NOSH,NOSH,NSXY,SXY)
end if
! NEXT, DETERMINE SEQUENTIAL SINGLE-TRANSFORMATION
! COEFFICIENTS TRA.
call PART(SXY,TRA1,TRA2)
call mma_deallocate(SXY)
if (PRTRA) then
  write(u6,*)
  call WRMAT('SINGLE-ORBITAL TRANSFORMATION COEFFS FOR STATE ONE (TRA1):',1,NOSH,NOSH,NTRA,TRA1)
  call WRMAT('FOR STATE TWO (TRA2):',1,NOSH,NOSH,NTRA,TRA2)
end if
! TRA1,TRA2 ARE NOT TRANSFORMATION MATRICES, ALTHOUGH THEY CONTAIN
! COEFFICIENTS FOR A SEQUENCE OF SINGLE-ORBITAL TRANSFORMATIONS.
! CONSTRUCT GENUINE TRANSFORMATION MATRICES (TEMPORARY) AND
! USE THEM TO TRANSFORM THE ORBITALS.
NCXA = NTRA
call mma_allocate(CXA,NCXA,Label='CXA')
call MKCXA(nIrrep,NOSH,NCXA,TRA1,CXA)
call TRAORB(nIrrep,NOSH,NBASF,NCXA,CXA,NCMO,CMO1)
call mma_deallocate(CXA)
NCYB = NTRA
call mma_allocate(CYB,NCYB,Label='CYB')
call MKCXA(nIrrep,NOSH,NCYB,TRA2,CYB)
call TRAORB(nIrrep,NOSH,NBASF,NCYB,CYB,NCMO,CMO2)
call mma_deallocate(CYB)
! print transformed MOs
if (PRORB) then
  write(u6,*)
  call WRMAT('TRANSFORMED MO COEFFICIENTS FOR STATE ONE (CMO1):',1,NBASF,NOSH,NCMO,CMO1)
  write(u6,*)
  call WRMAT('TRANSFORMED MO COEFFICIENTS FOR STATE TWO (CMO2):',1,NBASF,NOSH,NCMO,CMO2)
end if

end subroutine FINDT
