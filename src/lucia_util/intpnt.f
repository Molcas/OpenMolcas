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
* Copyright (C) 2000, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE INTPNT(IPNT1,ISL1,IPNT2,ISL2)
*
* Pointers to symmetry blocks of integrals
* IPNT1 : Pointer to given one-electron block, total symmetric
* ISL1  : Symmetry of last index for given first index, 1 e-
* IPNT2 : Pointer to given two-electron block
* ISL1  : Symmetry of last index for given first index, 1 e-
*
* In addition pointers to one-electron integrals with general
* symmetry is generated in WORK(KPGINT1(ISM))
*
* Pointers for similarity transformed Hamiltonian may also be
* generated
*
* Jeppe Olsen, Last Update : August 2000
      IMPLICIT REAL*8(A-H,O-Z)
*
* =====
*.Input
* =====
*
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "glbbas.fh"
#include "lucinp.fh"
#include "orbinp.fh"
#include "csm.fh"
#include "crun.fh"

*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
#include "cintfo.fh"
*
* =======
*. Output
* =======
*
      INTEGER IPNT1(NSMOB),ISL1(NSMOB)
      INTEGER IPNT2(NSMOB,NSMOB,NSMOB),ISL2(NSMOB,NSMOB,NSMOB)
*.0 : Pointers to one-integrals, all symmetries, Lower half matrices
      DO ISM = 1, NSMOB
        CALL PNT2DM(        1,    NSMOB,    NSMSX,    ADSXA,   NTOOBS,
     &                 NTOOBS,ISM  ,    ISL1,IWORK(KPGINT1(ISM)),MXPOBS)
      END DO
*.0.5 : Pointers to one-electron integrals, all symmetries, complete form
      DO ISM = 1, NSMOB
        CALL PNT2DM(        0,    NSMOB,    NSMSX,    ADSXA,   NTOOBS,
     &                 NTOOBS,ISM  ,   ISL1,IWORK(KPGINT1A(ISM)),MXPOBS)
      END DO
*.1 : Number of one-electron integrals
      CALL PNT2DM(        1,    NSMOB,    NSMSX,    ADSXA,   NTOOBS,
     &               NTOOBS,    ITSSX,     ISL1,    IPNT1,   MXPOBS)
*.2 : two-electron integrals
      CALL PNT4DM(    NSMOB,    NSMSX,   MXPOBS,   NTOOBS,   NTOOBS,
     &               NTOOBS,   NTOOBS,    ITSDX,    ADSXA,   SXDXSX,
     &                 I12S,     I34S,   I1234S,    IPNT2,     ISL2,
     &                ADASX)
*
c      IF(ISIMTRH.EQ.1) THEN
c*. Pointers for similarity transformed Hamiltonian
c        CALL PNT2DM(0,NSMOB,NSMSX,ADSXA,NTOOBS,NTOOBS,
c     &         1  ,ISL1,WORK(KPINT1_SIMTRH),MXPOBS)
c        I12 = 0
c        I34 = 0
c        I1234 = 1
c        CALL PNT4DM(NSMOB,NSMSX,MXPOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,
c     &              ITSDX,ADSXA,SXDXSX,I12,I34,I1234,
c     &              WORK(KPINT2_SIMTRH),ISL2,ADASX)
c      END IF
C?    write(6,*) ' Memory check INTPNT 2 '
C?    CALL MEMCHK
      RETURN
      END
************************************************************************
