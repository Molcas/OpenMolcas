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
      SUBROUTINE ZSYM1(NIRREP,IPRNT)
*
* Number of symmetries for d2h
* Symmetry connecting arrays
* ( trivial, written for compatibility with higher point groups)
*
      INTEGER SYMPRO(8,8)
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.fh"
*
C     PARAMETER ( MXPOBS = 20 )
#include "mxpdim.fh"
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

      NSMSX = NIRREP
      NSMDX = NIRREP
      NSMST = NIRREP
      NSMCI = NIRREP
      NSMXT = NIRREP
      ITSSX = 1
      ITSDX = 1
      ITSXT = 1

*
C     COPMT2(AIN,AOUT,NINR,NINC,NOUTR,NOUTC,IZERO)
      CALL ICPMT2(   SYMPRO,    ADASX,        8,        8,   MXPOBS,
     &               MXPOBS,        1)
      CALL ICPMT2(   SYMPRO,    ADSXA,        8,        8,   MXPOBS,
     &             2*MXPOBS,        1)
      CALL ICPMT2(   SYMPRO,    ASXAD,        8,        8,   MXPOBS,
     &             2*MXPOBS,        1)
      CALL ICPMT2(   SYMPRO,   SXSXDX,        8,        8, 2*MXPOBS,
     &             2*MXPOBS,        1)
      CALL ICPMT2(   SYMPRO,   SXDXSX,        8,        8, 2*MXPOBS,
     &             4*MXPOBS,        1)
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRNT)
      END
