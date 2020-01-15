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

      SUBROUTINE MKDYSZZ(CMOA,DYSAB,DYSZZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DYSAB(*),DYSZZ(*)
      DIMENSION CMOA(NCMO)
      INTEGER IBIO,IZZ,SYMOFF,BIOOFF
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
#include "rassi.fh"

C *** Re-express the DO coefficients in biorth basis DYSAB
C *** into atomic basis DYSZZ with help of CMOA that contains
C *** biorth orbitals in ZZ basis

      SYMOFF=0
      DO ISY1=1,NSYM
        NO1=NOSH(ISY1)
        NB1=NBASF(ISY1)
        DO IBIO=1,NO1
         DO IZZ=1,NB1
          BIOOFF=(IBIO-1)*NB1
          COEFF=DYSAB(IBIO)*CMOA(SYMOFF+BIOOFF+IZZ)
          DYSZZ(IZZ)=DYSZZ(IZZ)+COEFF
         END DO
        END DO
        SYMOFF=SYMOFF+NO1
      END DO

      RETURN
      END


