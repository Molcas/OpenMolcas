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
* Copyright (C) 2011, Jeppe Olsen                                      *
*               2011, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE ADDDIA_TERMS(NAEL,IASTR,NBEL,IBSTR,
     &                  NORB,CVEC,SVEC,NSMST,H,
     &                  XA,XB,SCR,RJ,RK,
     &                  NSSOA,NSSOB,
     &                  ECORE,
     &                  IPRNT,NTOOB,RJKAA,
     &                  IASPGP,IASM,IBSPGP,IBSM,FACTOR)
*
*. Update Sigma vector with diagonal terms for a given block
*     SVEC(IASPGP,IBSPGP) = SVEC(IASPGP,IBSPGP)
*                         + FACTOR*DIAG(IASPGP,IBSPGP)CVEC(IASPGP,IBSPGP)
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen and Giovanni Li Manni, September 2011
*
* I12 = 1 => only one-body part
*     = 2 =>      one+two-body part
      IMPLICIT REAL*8           (A-H,O-Z)
*. Input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION H(NORB)
      DIMENSION CVEC(*)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION SVEC(*)
*
      NTEST =  00
      NTEST = MAX(NTEST,IPRNT)
      IDUM = 0
      I12 = 2
c      IF(LUIN.GT.0) REWIND LUIN
c      IF(LUOUT.GT.0) REWIND LUOUT

      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' ADDDIA_TERMS in action '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' IASM, IASPGP, IBSM, IBSPGP = ',
     &               IASM, IASPGP, IBSM, IBSPGP
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMAT(H,1,NORB,1,NORB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
        END IF
        WRITE(6,*) ' FACTOR = ',FACTOR
      END IF
*
**3 Diagonal elements according to Handys formulae
*   (corrected for error)
*
*   DIAG(IDET) = HII*(NIA+NIB)
*              + 0.5 * ( J(I,J)-K(I,J) ) * N(I,A)*N(J,A)
*              + 0.5 * ( J(I,J)-K(I,J) ) * N(I,B)*N(J,B)
*              +         J(I,J) * N(I,A)*N(J,B)
* N(X) are occupation numbers
*
*. K goes to J - K
      IF(I12.EQ.2) CALL VECSUM(RK,RK,RJ,-1.0D0,+1.0D0,NTOOB **2)
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IASPGP
      CALL GETSTR_TOTSM_SPGP(1,IASPGP,IASM,NAEL,NASTR1,IASTR,
     &                         NORB,0,IDUM,IDUM)
*
      NIA = NSSOA(IASM,IASPGP)
*
      IF(NTEST.GE.1000) THEN
        write(6,*) ' After GETSTR for A strings '
        WRITE(6,*) ' alpha strings obtained '
        CALL IWRTMA(IASTR,NAEL,NIA,NAEL,NIA)
      END IF
*
      DO IA = 1 ,NIA
        EAA = 0.0D0
        DO IEL = 1, NAEL
          IAEL = IASTR(IEL,IA)
          EAA = EAA + H(IAEL)
          IF(I12.EQ.2) THEN
            DO JEL = 1, NAEL
              EAA =   EAA + 0.5D0*RK(IASTR(JEL,IA),IAEL )
            END DO
          END IF
        END DO
        RJKAA(IA) = EAA
      END DO
*. Obtain alpha strings of sym IBSM and type IBTP
      CALL GETSTR_TOTSM_SPGP(2,IBSPGP,IBSM,NBEL,NBSTR1,IBSTR,
     &                       NORB,0,IDUM,IDUM)
      NIB =  NSSOB(IBSM,IBSPGP)
      IDET = 0
      DO IB = 1 ,NIB
*. Terms depending only on IB
        HB = 0.0D0
        RJBB = 0.0D0
        CALL SETVEC(XB,0.0D0,NORB)

        DO IEL = 1, NBEL
          IBEL = IBSTR(IEL,IB)
          HB = HB + H(IBEL )
*
          IF(I12.EQ.2) THEN
            DO JEL = 1, NBEL
              RJBB = RJBB + RK(IBSTR(JEL,IB),IBEL )
            END DO
*
            DO IORB = 1, NORB
              XB(IORB) = XB(IORB) + RJ(IORB,IBEL)
            END DO
          END IF
        END DO

        EB = HB + 0.5D0*RJBB + ECORE
        DO IA = 1,NSSOA(IASM,IASPGP)
          IDET = IDET + 1
          X = EB + RJKAA(IA)
          DO IEL = 1, NAEL
            X = X +XB(IASTR(IEL,IA))
          END DO
          SVEC(IDET) = SVEC(IDET) + CVEC(IDET)*(X+FACTOR)
        END DO ! IA
      END DO ! IB
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input and output vectord, ADDDIA_TERMS '
        CALL WRTMAT(CVEC,1,IDET,1,IDET)
        CALL WRTMAT(SVEC,1,IDET,1,IDET)
      END IF
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_real_array(XA)
        CALL Unused_real_array(SCR)
      END IF
      END
