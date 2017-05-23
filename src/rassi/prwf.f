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
      SUBROUTINE PRWF(ISGSTRUCT,ICISTRUCT,ISYCI,CI,CITHR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CI(*)
#include "Struct.fh"
      Dimension iSGStruct(nSGSize)
      Dimension iCIStruct(nCISize)
#include "WrkSpc.fh"

      NSYM  =ISGSTRUCT(1)
      NLEV  =ISGSTRUCT(2)
      LISM  =ISGSTRUCT(3)
      NMIDV =ICISTRUCT(1)
      LNOW  =ICISTRUCT(3)
      LIOW  =ICISTRUCT(4)
      LNOCSF=ICISTRUCT(6)
      LIOCSF=ICISTRUCT(7)
      CALL GETMEM('ICS','ALLO','INTE',LICS,NLEV)
      CALL PRWF1(ISGSTRUCT,ICISTRUCT,NLEV,NMIDV,
     &           IWORK(LISM),IWORK(LICS),IWORK(LNOCSF),
     &           IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           ISYCI,CI,CITHR)
      CALL GETMEM('ICS','FREE','INTE',LICS,NLEV)
      RETURN
      END
