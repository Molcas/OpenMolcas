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
* Copyright (C) 2000, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE NONA2(
#define _CALLING_
#include "grd_mck_interface.fh"
     &               )
************************************************************************
* OBJECT: TO COMPUTE THE 2ND DERIVATIVE NONADIABATIC COUPLING
* INTEGRALS, OF TYPE
*     < D/DX CHI_1 | D/DX CHI_2 >
*
*     AUTHOR: PER AKE MALMQVIST, MAX PLANCK INSTITUT F ASTROPHYSIK
*             GARCHING, MUENCHEN NOV 2000
*     AFTER PROGRAMMING PATTERN ESTABLISHED BY ROLAND LINDH
*
************************************************************************
      use Her_RW
      use Center_Info
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"

#include "grd_mck_interface.fh"

*     Local variables

      LOGICAL ABEQ(3)
C The following call parameters are not used:
C IDCNT,ISTABM,NSTABM,ZINV
C They must still be present, because the call parameter list must

      NELEM(LA)=(LA+2)*(LA+1)/2
      ABEQ(1) = A(1).EQ.RB(1)
      ABEQ(2) = A(2).EQ.RB(2)
      ABEQ(3) = A(3).EQ.RB(3)

      NIP = 1
      IPAXYZ = NIP
      NIP = NIP + NZETA*3*NHER*(LA+2)
      IPBXYZ = NIP
      NIP = NIP + NZETA*3*NHER*(LB+2)
      IPRXYZ = NIP
      NIP = NIP + NZETA*3*NHER*(NORDOP+1)
      IPRNXYZ = NIP
      NIP = NIP + NZETA*3*(LA+2)*(LB+2)*(NORDOP+1)
      IPALPH = NIP
      NIP = NIP + NZETA
      IPBETA = NIP
      NIP = NIP + NZETA
      IPSCRT=NIP
      NIP=NIP+NELEM(LA)*NELEM(LB)*NZETA*2


      IF (NIP-1.GT.NARR) THEN
        WRITE(6,*)' NONA2: Too small array.'
        WRITE(6,*)' Submitted array size NARR=',NARR
        WRITE(6,*)' Needed size at least NIP =',NIP
        CALL Abend
      END IF

* COMPUTE THE CARTESIAN VALUES OF THE BASIS FUNCTIONS ANGULAR PART
      CALL CRTCMP(ZETA,P,NZETA,A,ARRAY(IPAXYZ),
     &               LA+1,HerR(iHerR(NHER)),NHER,ABEQ)
      CALL CRTCMP(ZETA,P,NZETA,RB,ARRAY(IPBXYZ),
     &               LB+1,HerR(iHerR(NHER)),NHER,ABEQ)

CPAM: WILL WE NEED THIS??
* COMPUTE THE CONTRIBUTION FROM THE MULTIPOLE MOMENT OPERATOR
      ABEQ(1) = .FALSE.
      ABEQ(2) = .FALSE.
      ABEQ(3) = .FALSE.
      CALL CRTCMP(ZETA,P,NZETA,CCOOR,ARRAY(IPRXYZ),
     &            NORDOP,HerR(iHerR(NHER)),NHER,ABEQ)

* COMPUTE THE PRIMITIVE 1-DIMENSIONAL OVERLAP INTEGRALS.
       CALL ASSMBL(ARRAY(IPRNXYZ),
     &             ARRAY(IPAXYZ),LA+1,
     &             ARRAY(IPRXYZ),NORDOP,
     &             ARRAY(IPBXYZ),LB+1,
     &             NZETA,HerR(iHerW(NHER)),NHER)

* COMBINE THE CARTESIAN COMPONENTS OF THE 2DC MATRIX ELEMENTS
      IP = IPALPH
      DO IBETA = 1, NBETA
         CALL DCOPY_(NALPHA,ALPHA,1,ARRAY(IP),1)
         IP = IP + NALPHA
      END DO
      IP = IPBETA
      DO IALPHA = 1, NALPHA
         CALL DCOPY_(NBETA,BETA,1,ARRAY(IP),NALPHA)
         IP = IP + 1
      END DO
      CALL CMBN2DC(ARRAY(IPRNXYZ),NZETA,LA,LB,ZETA,
     &            RKAPPA,ARRAY(IPSCRT),
     &            ARRAY(IPALPH),ARRAY(IPBETA),
     &            IFGRAD)

* SYMMETRY ADAPT THE 2ND DERIVATIVE COUPLING INTEGRALS
      CALL SYMADO_MCK(ARRAY(IPSCRT),NZETA*NELEM(LA)*NELEM(LB),
     &            FINAL,NROP,
     &            nOP,LOPER,INDGRD,IU,IV,IFGRAD,IDCAR,TRANS)

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_real_array(ZINV)
         CALL Unused_integer(IDCNT)
         CALL Unused_integer_array(ISTABM)
         CALL Unused_integer(NSTABM)
      END IF
      END
