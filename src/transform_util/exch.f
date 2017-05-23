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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE EXCH(ISYP,ISYI,ISYQ,ISYJ,II,IJ,ERI,SCR)

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION ERI(*), SCR(*)
      LOGICAL TRANSP
#include "intgrl.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "blksize.fh"

      NDIM2M=(NSYMZ*(NSYMZ+1))/2
      IF(ISYP.GE.ISYQ) THEN
        ISY1=ISYP
        ISY2=ISYQ
        ISY3=ISYI
        ISY4=ISYJ
        I3=II
        I4=IJ
        TRANSP=.TRUE.
      ELSE
        ISY1=ISYQ
        ISY2=ISYP
        ISY3=ISYJ
        ISY4=ISYI
        I3=IJ
        I4=II
        TRANSP=.FALSE.
      END IF

C We are going for an EXCH1 or EXCH2 integral of symmetry
C block ISY1,ISY3,ISY2,ISY4.
C Disk address to symmetry block, and record nr:

      IS12=(ISY1*(ISY1-1))/2+ISY2
      IF(ISY3.GT.ISY4) THEN
        IS34=(ISY3*(ISY3-1))/2+ISY4
        IDISK=IAD2M(2,IS12+NDIM2M*(IS34-1))
        IREC=I4+NOSHZ(ISY4)*(I3-1)
      ELSE IF(ISY3.EQ.ISY4) THEN
        IS34=(ISY3*(ISY3-1))/2+ISY4
        IF(I3.GE.I4) THEN
          IDISK=IAD2M(2,IS12+NDIM2M*(IS34-1))
          IREC=(I3*(I3-1))/2 + I4
        ELSE
          IDISK=IAD2M(2,IS12+NDIM2M*(IS34-1))
          IREC=(I4*(I4-1))/2 + I3
          TRANSP= .NOT. TRANSP
        END IF
      ELSE
        IS34=(ISY4*(ISY4-1))/2+ISY3
        IDISK=IAD2M(3,IS12+NDIM2M*(IS34-1))
        IREC=I3+NOSHZ(ISY3)*(I4-1)
      END IF

C Buffer size:
      NO1=NORBZ(ISY1)
      NO2=NORBZ(ISY2)
      NBUF=NO1*NO2
      IF(NBUF.EQ.0) RETURN

C Address update for earlier records, then read:
* PAM07 * Eliminate unsafe IPOSFILE call
*      IDISK=IDISK+iPosFile(NBUF)*(IREC-1)
* Replace with dummy i/o operations: Is this efficience issue??
      DO I=1,IREC-1
        CALL dDAFILE(LUINTMZ,0,SCR,NBUF,IDISK)
      END DO

      IF(TRANSP) THEN
        CALL dDAFILE(LUINTMZ,2,SCR,NBUF,IDISK)
        CALL TRNSPS(NO2,NO1,SCR,ERI)
      ELSE
        CALL dDAFILE(LUINTMZ,2,ERI,NBUF,IDISK)
      END IF

      RETURN
      END
