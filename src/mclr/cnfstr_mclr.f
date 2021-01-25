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
* Copyright (C) 1989,1991,1993, Jeppe Olsen                            *
************************************************************************
      SUBROUTINE CNFSTR_MCLR(ICONF,ITYP,IASTR,IBSTR,NORB,NAEL,NBEL,
     &                  IDET,IPRODT,IAGRP,IBGRP,ISCR,SIGN,NTEST)
*
* An orbital configuration ICONF is given,
* Obtain the corresponding alpha strings,IASTR
*        the corresponding beta  strings,IBSTR
*        the corresponding sign array   ,ISIGN
*
* Jeppe Olsen , Summer of 89
*
* Ninob added , October 1991
*
* Modified September 1993 for LUCIA
*
      IMPLICIT REAL*8 ( A-H,O-Z)
*. Specific input
      INTEGER ICONF(*)
*. General input
      INTEGER IPRODT(*)
*
*. Scratch : required length : IDET * ( NAEL + NBEL ) + 2(NAEL+NBEL)
* ( this includes NAEL+NBEL words needed in DETSTR )
      INTEGER ISCR(*)
*. Output
      DIMENSION IASTR(*),IBSTR(*), SIGN(*)
*
#include "detdim.fh"
#include "spinfo_mclr.fh"
*
*
      NEL = NAEL + NBEL
      IOPEN = ITYP-1+MINOP
      ICLOS = (NAEL + NBEL - IOPEN) / 2
*
** Spin orbital occupations of determinants of configuration
*
      KLFREE = 1
*
      KLDETS = KLFREE
      KLFREE = KLDETS + IDET*(NAEL+NBEL)
*
      KLIA = KLFREE
      KLFREE = KLFREE + NAEL
*
      KLIB = KLFREE
      KLFREE = KLFREE + NBEL
*
      KLDET  = KLFREE
      KLFREE = KLFREE + NAEL + NBEL
*

* Pointer for determinants of prototype ITYP

      IP = 1
      DO 10 JTYP = 1, ITYP - 1
        IP = IP + NDPCNT(JTYP)*(JTYP-1+MINOP)
   10 CONTINUE
*. Expand into determinants
*
      CALL CNDET_MCLR(ICONF,IPRODT(IP),IDET,NAEL+NBEL,NORB,
     &           IOPEN,ICLOS,ISCR(KLDETS),NTEST)
*
** Separate determinants into strings and determine sign change
*
      DO 100 JDET = 1, IDET
        CALL DETSTR_MCLR(ISCR(KLDETS+(JDET-1)*NEL),
     &              ISCR(KLIA),ISCR(KLIB),
     &              NEL,NAEL,NBEL,
     &              NORB,ISIGN,ISCR(KLDET),NTEST)
*. Actual numbers of alpha and beta string
        IASTR(JDET) = ISTRN_MCLR(ISCR(KLIA),IAGRP)
        IBSTR(JDET) = ISTRN_MCLR(ISCR(KLIB),IBGRP)
        SIGN(JDET) = DBLE(ISIGN)
  100 CONTINUE
*
*
      RETURN
      END
