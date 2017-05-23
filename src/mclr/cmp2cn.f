************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE CMP2CN(ICNL,NCLL,NOPL,ICNR,NCLR,NOPR,ISCR,NORB,NDIFF,
     &                  NTEST)
*
* Number of differences in occupation of two configurations
*
      DIMENSION ICNL(*),ICNR(*)
      DIMENSION ISCR(*)
* Length of Scratch : Number of orbitals
*
      CALL ISETVC(ISCR,0,NORB)
      DO 10 ICL = 1, NCLL
        ISCR(ICNL(ICL) ) = 2
   10 CONTINUE
      DO 20 IOP = 1, NOPL
        ISCR(ICNL(NCLL+IOP)) = 1
   20 CONTINUE
*
      NDIFF = 0
      DO 30 ICL = 1, NCLR
        NDIFF = NDIFF +  2 - ISCR(ICNR(ICL))
   30 CONTINUE
      DO 40 IOP = 1, NOPR
        IF(ISCR(ICNR(NCLR+IOP)).EQ.0) NDIFF = NDIFF + 1
   40 CONTINUE
*
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NTEST)
      END
