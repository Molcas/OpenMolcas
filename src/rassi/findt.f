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
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CMO1(NCMO),CMO2(NCMO)
      DIMENSION TRA1(NTRA),TRA2(NTRA)
#include "rasdim.fh"
#include "cntrl.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "WrkSpc.fh"
C BEGIN BY COMPUTING MO OVERLAPS, SXY.
C MO OVERLAP MATRIX HAS SAME STRUCTURE AS TRA1,TRA2:
      NSXY=NTRA
      CALL GETMEM('SXY   ','ALLO','REAL',LSXY,NSXY)
      CALL MKSXY(CMO1,CMO2,WORK(LSXY))
      IF (PRSXY) THEN
        WRITE(6,*)
        CALL WRMAT('MO OVERLAP MATRIX IN ORIGINAL ORBITAL BASES:',
     *             1,NOSH,NOSH,NSXY,WORK(LSXY))
      END IF
C NEXT, DETERMINE SEQUENTIAL SINGLE-TRANSFORMATION
C COEFFICIENTS TRA.
      CALL PART(WORK(LSXY),TRA1,TRA2)
      CALL GETMEM('      ','FREE','REAL',LSXY,NSXY)
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
      CALL GETMEM('CXA   ','ALLO','REAL',LCXA,NCXA)
C CALL MKCXA(NSYM,NOSH,NCXA,TRA,CXA)
      CALL MKCXA(NSYM,NOSH,NCXA,TRA1,WORK(LCXA))
C CALL TRAORB(NSYM,NOSH,NBASF,NCXA,CXA,NCMO,CMO)
      CALL TRAORB(NSYM,NOSH,NBASF,NCXA,WORK(LCXA),NCMO,CMO1)
      CALL GETMEM('      ','FREE','REAL',LCXA,NCXA)
      NCYB=NTRA
      CALL GETMEM('CYB   ','ALLO','REAL',LCYB,NCYB)
      CALL MKCXA(NSYM,NOSH,NCYB,TRA2,WORK(LCYB))
      CALL TRAORB(NSYM,NOSH,NBASF,NCYB,WORK(LCYB),NCMO,CMO2)
      CALL GETMEM('      ','FREE','REAL',LCYB,NCYB)
      RETURN
      END
