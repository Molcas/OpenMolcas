************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1989, Per Ake Malmqvist                                *
************************************************************************
*****************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST
*  SUBROUTINE FINDT     IBM-3090 RELEASE 89 01 31
*  GIVEN CMO1 AND CMO2, FIND SUITABLE COUEFFICIENTS FOR THE
*  SINGLE-ORBITAL TRANSFORMATION SEQUENCE FOR EACH SET, RETURNED
*  IN ARRAYS TRA1 AND TRA2, FOR TRANSFORMATION INTO BIORTHONORMAL
*  ORBITALS.
*****************************************************************
      SUBROUTINE FINDT (CMO1,CMO2,TRA1,TRA2)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: PRSXY, PRTRA, PRORB
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: NCMO,NTRA,NSXY,NBASF,NOSH
      IMPLICIT None
      Real*8 CMO1(NCMO),CMO2(NCMO)
      Real*8 TRA1(NTRA),TRA2(NTRA)
      Real*8, allocatable:: SXY(:), CXA(:), CYB(:)
      Integer NCXA, NCYB

C BEGIN BY COMPUTING MO OVERLAPS, SXY.
C MO OVERLAP MATRIX HAS SAME STRUCTURE AS TRA1,TRA2:
      NSXY=NTRA
      CALL mma_allocate(SXY,NSXY,Label='SXY')
      CALL MKSXY(CMO1,CMO2,SXY)
      IF (PRSXY) THEN
        WRITE(6,*)
        CALL WRMAT('MO OVERLAP MATRIX IN ORIGINAL ORBITAL BASES:',
     *             1,NOSH,NOSH,NSXY,SXY)
      END IF
C NEXT, DETERMINE SEQUENTIAL SINGLE-TRANSFORMATION
C COEFFICIENTS TRA.
      CALL PART(SXY,TRA1,TRA2)
      CALL mma_deallocate(SXY)
      IF (PRTRA) THEN
        WRITE(6,*)
        CALL WRMAT(
     *  'SINGLE-ORBITAL TRANSFORMATION COEFFS FOR STATE ONE (TRA1):',
     *             1,NOSH,NOSH,NTRA,TRA1)
        CALL WRMAT('FOR STATE TWO (TRA2):',
     *             1,NOSH,NOSH,NTRA,TRA2)
      END IF
C TRA1,TRA2 ARE NOT TRANSFORMATION MATRICES, ALTHOUGH THEY CONTAIN
C COEFFICIENTS FOR A SEQUENCE OF SINGLE-ORBITAL TRANSFORMATIONS.
C CONSTRUCT GENUINE TRANSFORMATION MATRICES (TEMPORARY) AND
C USE THEM TO TRANSFORM THE ORBITALS.
      NCXA=NTRA
      CALL mma_allocate(CXA,NCXA,Label='CXA')
      CALL MKCXA(NSYM,NOSH,NCXA,TRA1,CXA)
      CALL TRAORB(NSYM,NOSH,NBASF,NCXA,CXA,NCMO,CMO1)
      CALL mma_deallocate(CXA)
      NCYB=NTRA
      CALL mma_allocate(CYB,NCYB,Label='CYB')
      CALL MKCXA(NSYM,NOSH,NCYB,TRA2,CYB)
      CALL TRAORB(NSYM,NOSH,NBASF,NCYB,CYB,NCMO,CMO2)
      CALL mma_deallocate(CYB)
C print transformed MOs
      if (PRORB) then
        write(6,*)
        call WRMAT('TRANSFORMED MO COEFFICIENTS FOR STATE ONE (CMO1):',
     &               1,NBASF,NOSH,NCMO,CMO1)
        write(6,*)
        call WRMAT('TRANSFORMED MO COEFFICIENTS FOR STATE TWO (CMO2):',
     &               1,NBASF,NOSH,NCMO,CMO2)
      end if

      END SUBROUTINE FINDT
