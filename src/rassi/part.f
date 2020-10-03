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
* Copyright (C) 1984,1989, Per Ake Malmqvist                           *
************************************************************************
      SUBROUTINE PART(SXY,TRA1,TRA2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SXY(NSXY),TRA1(NTRA),TRA2(NTRA)
      DIMENSION NSIZE(4)
C  PURPOSE: SXY CONTAINS THE NONSECONDARY PART OF THE MO OVERLAP
C  MATRIX. UPON RETURN, TRA1 AND TRA2 WILL CONTAIN THE COEFFICIENTS
C  FOR SEQUENTIAL SINGLE-ORBITAL TRANSFORMATIONS (VI.2, MY IJQC ARTICLE)
C  TO BIORTHONORMAL ORBITALS. SXY, TRA1 AND TRA2 ARE SYMMETRY-BLOCKED.
C  ORIGINAL VERSION, MALMQUIST 84-04-04
C  RASSCF VERSION,   MALMQUIST 89-11-15
#include "WrkSpc.fh"
#include "rasdef.fh"
#include "symmul.fh"
#include "rassi.fh"
      NOMAX=0
      DO 5 ISY=1,NSYM
        NOMAX=MAX(NOSH(ISY),NOMAX)
5     CONTINUE
      CALL GETMEM('SCRMAT','ALLO','REAL',LSCRMAT,NOMAX*NOMAX)
      CALL GETMEM('SCRPIV','ALLO','INTE',LSCRPIV,2*NOMAX)
      CALL GETMEM('SCRBUF','ALLO','REAL',LSCRBUF,NOMAX)
      II=1
      DO 10 ISY=1,NSYM
        NDIMEN=NOSH(ISY)
        IF(NDIMEN.EQ.0) GOTO 10
        NBLOCK=0
        N=NISH(ISY)
        IF(N.GT.0) THEN
          NBLOCK=NBLOCK+1
          NSIZE(NBLOCK)=N
        END IF
        N=NRS1(ISY)
        IF(N.GT.0) THEN
          NBLOCK=NBLOCK+1
          NSIZE(NBLOCK)=N
        END IF
        N=NRS2(ISY)
        IF(N.GT.0) THEN
          NBLOCK=NBLOCK+1
          NSIZE(NBLOCK)=N
        END IF
        N=NRS3(ISY)
        IF(N.GT.0) THEN
          NBLOCK=NBLOCK+1
          NSIZE(NBLOCK)=N
        END IF
        CALL PART1(NDIMEN,NBLOCK,NSIZE,SXY(II),TRA1(II),TRA2(II),
     *             WORK(LSCRMAT),IWORK(LSCRPIV),WORK(LSCRBUF))
        II=II+NDIMEN**2
10      CONTINUE
      CALL GETMEM('SCRMAT','FREE','REAL',LSCRMAT,NOMAX*NOMAX)
      CALL GETMEM('SCRPIV','FREE','INTE',LSCRPIV,2*NOMAX)
      CALL GETMEM('SCRBUF','FREE','REAL',LSCRBUF,NOMAX)
      RETURN
      END
