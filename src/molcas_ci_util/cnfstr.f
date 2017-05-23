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
* Copyright (C) 1989, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE CNFSTR(ICONF,ITYP,IASTR,IBSTR,NORB,NAEL,NBEL,
     *                  IDET,IPRODT,ISCR, SIGN,IPREXH)
*
* An orbital configuration ICONF is given,
* Obtain the corresponding alpha strings,IASTR
*        the corresponding beta  strings,IBSTR
*        the corresponding sign array   ,ISIGN
*
* Jeppe Olsen , Summer of '89
*
      IMPLICIT REAL*8 ( A-H,O-Z)
      DIMENSION ICONF(*),ISCR(*),IPRODT(*)
*PAM06 NOTE: NAEL and NBEL can legally be=0.
*PAM06      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*), SIGN(*)
      DIMENSION IASTR(*),IBSTR(*), SIGN(*)
#include "spinfo.fh"
#include "ciinfo.fh"
*
      NEL = NAEL + NBEL
      IOPEN = ITYP-1+MINOP
      ICLOS = (NAEL + NBEL - IOPEN) / 2
      IOCC  = IOPEN + ICLOS
*
** Spin orbital occupations of determinants of configuration
*
      KLFREE = 1
      KLDETS = KLFREE
      KLFREE = KLDETS + IDET*(NAEL+NBEL)
* Pointer for determinants of prototype ITYP
      IP = 1
      DO 10 JTYP = 1, ITYP - 1
        IP = IP + NDTFTP(JTYP)*(JTYP-1+MINOP)
10    CONTINUE
      CALL CNDET(ICONF,IPRODT(IP),IDET,NAEL+NBEL,NORB,
     &           IOPEN,ICLOS,ISCR(KLDETS),IPREXH)
*
** Separate determinants into strings and determine sign change
*
      DO 100 JDET = 1, IDET
        CALL DETSTR2(ISCR(KLDETS+(JDET-1)*NEL),
*PAM06     &              IASTR(1,JDET),IBSTR(1,JDET),
     &        IASTR(1+NAEL*(JDET-1)),IBSTR(1+NBEL*(JDET-1)),
     &              NEL,NAEL,NBEL,ISIGN,ISCR(KLFREE),IPREXH)
        SIGN(JDET) = DBLE(ISIGN)
100   CONTINUE
*
      NTEST =  0
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Output from CNFSTR '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' Input configuration '
        CALL IWRTMA(ICONF,1,IOCC,1,IOCC)
        WRITE(6,*) ' Corresponding alpha and beta strings'
        CALL IWRTMA(IASTR,NAEL,IDET,NAEL,IDET)
        CALL IWRTMA(IBSTR,NBEL,IDET,NBEL,IDET)
        WRITE(6,*) ' SIGN ARRAY '
        CALL WRTMAT( SIGN,1,IDET,1,IDET)
      END IF
*
      RETURN
      END
