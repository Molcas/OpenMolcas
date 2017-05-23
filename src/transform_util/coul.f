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
      SUBROUTINE COUL(ISYP,ISYQ,ISYI,ISYJ,II,IJ,ERI,SCR)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "standard_iounits.fh"
      DIMENSION ERI(*), SCR(*)
      LOGICAL TRANSP
      LOGICAL TRIANG

#include "intgrl.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "blksize.fh"

C Return a matrix ERI(IP,IQ) of two-electron integrals (pq,ij).
C IP=1..NORB(ISYP), IQ=1..NORB(ISYQ),1<=II<=NOSH(ISYI),
C 1<=IJ<=NOSH(ISYJ). Normal Fortran storage order.

      NDIM2M=(NSYMZ*(NSYMZ+1))/2
        TRANSP=.false.
        TRIANG=.false.
      IF(ISYP.GE.ISYQ) THEN
        ISY1=ISYP
        ISY2=ISYQ
        TRANSP=.true.
c        TYPE='TRANSPOS'
        IF(ISYP.EQ.ISYQ) then
        TRIANG=.true.
        endif
c        IF(ISYP.EQ.ISYQ) TYPE='TRIANGUL'
      ELSE
        ISY1=ISYQ
        ISY2=ISYP
c        TYPE='NORMAL  '
      END IF
      IF(ISYI.GE.ISYJ) THEN
        ISY3=ISYI
        ISY4=ISYJ
        I3=II
        I4=IJ
        IF(ISYI.EQ.ISYJ .AND. II.LT.IJ) THEN
          I3=IJ
          I4=II
        END IF
      ELSE
        ISY3=ISYJ
        ISY4=ISYI
        I3=IJ
        I4=II
      END IF

C Disk address to symmetry block:
      IS12=(ISY1*(ISY1-1))/2+ISY2
      IS34=(ISY3*(ISY3-1))/2+ISY4
      IDISK=IAD2M(1,IS12+NDIM2M*(IS34-1))

C Record number, in symmetry block:
      I34=I4+NOSHZ(ISY4)*(I3-1)
      IF(ISY3.EQ.ISY4) I34=(I3*(I3-1))/2+I4

C Buffer sizes:
      NO1=NORBZ(ISY1)
      NO2=NORBZ(ISY2)
      NBUF=NO1*NO2
      IF(TRIANG) NBUF=(NBUF+NO1)/2
      IF(NBUF.EQ.0) RETURN

C Address update for earlier records, then read:
      IDISK=IDISK+NBUF*(I34-1)
*
      IF(.not.TRANSP) THEN
        CALL dDAFILE(LUINTMZ,2,ERI,NBUF,IDISK)
      ELSE
        CALL dDAFILE(LUINTMZ,2,SCR,NBUF,IDISK)
*
        IF(.not.TRIANG) THEN
          CALL TRNSPS(NO2,NO1,SCR,ERI)
        ELSE
          CALL SQUARE(SCR,ERI,NO1,1,NO1)
        END IF
      END IF
*
      RETURN
      END
