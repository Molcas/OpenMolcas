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
* Copyright (C) 2018, Jesper Norell                                    *
************************************************************************

*****************************************************************
*  SUBROUTINE MKDYSZZ
*  PURPOSE: CALCULATE DYSON ORBITAL COEFFICIENTS FOR CI EXPANSIONS IN
*  BASIS FUNCTION BASE BASE Z,
*  IN ANALOGUE TO MKTDZZ FOR TRANSITION DENSITY MATRIX.
*****************************************************************
*  MODIFIED BY BRUNO TENORIO TO ADDRESS SYMMETRY
*  SEPTEMBER 2020
*****************************************************************
      SUBROUTINE MKDYSZZ(CMOA,DYSAB,DYSZZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DYSAB(*),DYSZZ(*)
      DIMENSION CMOA(NCMO)
      INTEGER IBIO,IZZ,SYMOFF,BIOOFF,IBIOFF
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
#include "rassi.fh"

C *** Re-express the DO coefficients in biorth basis DYSAB
C *** into atomic basis DYSZZ with help of CMOA that contains
C *** biorth orbitals in ZZ basis

      SYMOFF=0
      IBIOFF=0
      IZZOFF=0
      DO ISY1=1,NSYM
        NO1=NOSH(ISY1)
        NA1=NASH(ISY1)
        NB1=NBASF(ISY1)
        IF(NA1.GT.0) THEN
        DO IBIO=1,NO1
         DO IZZ=1,NB1
          BIOOFF=(IBIO-1)*NB1
          COEFF=DYSAB(IBIO+IBIOFF)*CMOA(SYMOFF+BIOOFF+IZZ)
          DYSZZ(IZZ+IZZOFF)=DYSZZ(IZZ+IZZOFF)+COEFF
         END DO
        END DO
        END IF
        IZZOFF=NB1+IZZOFF
        IBIOFF=NO1+IBIOFF
        SYMOFF=(NO1*NB1)+SYMOFF
      END DO


      RETURN
      END
