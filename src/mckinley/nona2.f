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
C The call from DRVH1_MCK is:
C            Call Cnt1El(NONA2,NA2Mem,Label,idcnt,idcar,loper,
C     &               One,.false.,Work(ipFock),
C     &               'NONA2   ',0)
C It is received by DRVH1_MCK as:
C      SubRoutine Cnt1El(Kernel,KrnlMm,Label,
C     &                 iDCnt,iDCar,rHrmt,DiffOp,dens,
C     &                 Lab_Dsk,iadd)
C The call is then transferred to NONA2 as:
C             Call Kernel(Work(iExp),iPrim,Work(jExp),jPrim,
C     &                   Work(iZeta),Work(ipZI),
C     &                   Work(iKappa),Work(iPCoor),
C     &                   Work(ipFnl),iPrim*jPrim,
C     &                   iAng,jAng,A,RB,nOrder,Work(iKern),
C     &                   MemKrn,Ccoor,nOrdOp,IfGrd,IndGrd,nop,
C     &                   loper,nStab(mdc+iCnt),
C     &                   nStab(ndc+jCnt),nic,idcar,idcnt,
C     &                   iStabM,nStabM,trans)
      SUBROUTINE NONA2(ALPHA,NALPHA,BETA, NBETA,ZETA,ZINV,RKAPPA,
     &               PCENT,FINAL,NZETA,LA,LB,ACENT,BCENT,NHER,
     &               ARRAY,NARR,CCOOR,NORDOP,IFGRAD,
     &               INDGRD,nOP,
     &               LOPER,IU,IV,NROP,IDCAR,IDCNT,ISTABM,NSTABM,TRANS)
************************************************************************
* OBJECT: TO COMPUTE THE 2ND DERIVATIVE NONADIABATIC COUPLING
* INTEGRALS, OF TYPE
*     < D/DX CHI_1 | D/DX CHI_2 >
*
* CALLED FROM: ONEEL
*
* CALLING    : QENTER
*              RECPRT
*              CRTCMP
*              ASSMBL
*              GETMEM
*              DCOPY   (ESSL)
*              CMBN2DC
*              QEXIT
*
*     AUTHOR: PER AKE MALMQVIST, MAX PLANCK INSTITUT F ASTROPHYSIK
*             GARCHING, MUENCHEN NOV 2000
*     AFTER PROGRAMMING PATTERN ESTABLISHED BY ROLAND LINDH
*
************************************************************************
      use Her_RW
C conform to that of all 'KERNEL' routines.
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
      INTEGER INDGRD(NIRREP), nOP(2)
      REAL*8 FINAL(NZETA,(LA+1)*(LA+2)/2,(LB+1)*(LB+2)/2,NROP),
     &       ZETA(NZETA), ZINV(NZETA), ALPHA(NALPHA), BETA(NBETA),
     &       RKAPPA(NZETA),PCENT(NZETA,3), ACENT(3), BCENT(3),CCOOR(3),
     &       ARRAY(NARR)
      LOGICAL ABEQ(3), IFGRAD(3,2),TRANS(2)
C The following call parameters are not used:
C IDCNT,ISTABM,NSTABM,ZINV
C They must still be present, because the call parameter list must

      NELEM(LA)=(LA+2)*(LA+1)/2
      ABEQ(1) = ACENT(1).EQ.BCENT(1)
      ABEQ(2) = ACENT(2).EQ.BCENT(2)
      ABEQ(3) = ACENT(3).EQ.BCENT(3)

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
      CALL CRTCMP(ZETA,PCENT,NZETA,ACENT,ARRAY(IPAXYZ),
     &               LA+1,HerR(iHerR(NHER)),NHER,ABEQ)
      CALL CRTCMP(ZETA,PCENT,NZETA,BCENT,ARRAY(IPBXYZ),
     &               LB+1,HerR(iHerR(NHER)),NHER,ABEQ)

CPAM: WILL WE NEED THIS??
* COMPUTE THE CONTRIBUTION FROM THE MULTIPOLE MOMENT OPERATOR
      ABEQ(1) = .FALSE.
      ABEQ(2) = .FALSE.
      ABEQ(3) = .FALSE.
      CALL CRTCMP(ZETA,PCENT,NZETA,CCOOR,ARRAY(IPRXYZ),
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
         CALL Unused_integer(ISTABM)
         CALL Unused_integer(NSTABM)
      END IF
      END
