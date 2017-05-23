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
* Copyright (C) 1991,1994, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE CIDIA4(NAEL,IASTR,NBEL,IBSTR,
     &                  NORB,DIAG,NSMST,H,
     &                  ISMOST,IBLTP,XA,XB,SCR,RJ,RK,
     &                  NSSOA,NSSOB,IOCOC,NOCTPA,NOCTPB,
     &                  ISSOA,ISSOB,LUDIA,ECORE,
     &                  PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR)
*
* Calculate determinant diagonal
* Turbo-ras version
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, Winter of 1991
* K => J - K moved outside, April 1994
*
      IMPLICIT REAL*8 (A-H,O-Z)
*.General input
      DIMENSION NSSOA(NOCTPA,*),NSSOB(NOCTPB,* )
      DIMENSION ISSOA(NOCTPA,*),ISSOB(NOCTPB,*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION H(NORB)
*. Specific input
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ISMOST(*),IBLTP(*)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
*. Output
      DIMENSION DIAG(*)
*
*
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0d0
      ELSE
         XADD = 0.0D0
      END IF
*
*
*
**3 Diagonal elements according to Handys formulae
*   (corrected for error)
*
*   DIAG(IDET) = HII*(NIA+NIB)
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIA*NJA
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIB*NJB
*              +         J(I,J) * NIA*NJB
*
*. K goes to J - K
*
      IDET = 0
      ITDET = 0
      IF(LUDIA.NE.0) REWIND LUDIA
      DO 1000 IASM = 1, NSMST
        IBSM = ISMOST(IASM)
        IF(IBSM.EQ.0.OR.IBLTP(IASM).EQ.0) GOTO 1000
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*
        DO 999  IATP = 1,NOCTPA
          IF(IREST1.EQ.1) THEN
            MXBTP = IATP
          ELSE
            MXBTP = NOCTPB
          END IF
          DO 900 IBTP = 1,MXBTP
          IF(IOCOC(IATP,IBTP) .EQ. 0 ) GOTO 900
          IBSTRT = ISSOB(IBTP,IBSM)
          IBSTOP = IBSTRT + NSSOB(IBTP,IBSM)-1
          DO 899 IB = IBSTRT,IBSTOP
            IBREL = IB - IBSTRT + 1
*
*. Terms depending only on IB
*
            ZERO = 0.0D0
*           CALL SETVEC(XB,ZERO,NORB)
            call dcopy_(NORB,ZERO,0,XB,1)
            HB = 0.0D0
            RJBB = 0.0D0
*           CALL SETVEC(XB,ZERO,NORB)
*
            DO 990 IEL = 1, NBEL
              IBEL = IBSTR(IEL,IB)
              HB = HB + H(IBEL )
*
              DO 980 JEL = 1, NBEL
                RJBB = RJBB + RK(IBSTR(JEL,IB),IBEL )
  980         CONTINUE
*
              DO 970 IORB = 1, NORB
                XB(IORB) = XB(IORB) + RJ(IORB,IBEL)
  970         CONTINUE
  990       CONTINUE
            EB = HB + 0.5D0*RJBB + ECORE
*
            IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
              IASTRT = ISSOA(IATP,IASM) - 1 + IBREL
            ELSE
              IASTRT = ISSOA(IATP,IASM)
            END IF
            IASTOP = ISSOA(IATP,IASM) + NSSOA(IATP,IASM) - 1
            DO 800 IA = IASTRT,IASTOP
              IDET = IDET + 1
              ITDET = ITDET + 1
              X1 = EB
              X2 = 0.0D0
              DO 890 IEL = 1, NAEL
                IAEL = IASTR(IEL,IA)
                X1 = X1 + ( H(IAEL )+XB(IAEL) )
                DO 880 JEL = 1, NAEL
                 X2 = X2 + RK(IASTR(JEL,IA),IAEL )
  880           CONTINUE
  890         CONTINUE
              DIAG(IDET) = X1 + 0.5D0*X2
              IF(IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
  800       CONTINUE
  899     CONTINUE
*. Yet a RAS block of the diagonal has been constructed
          IF(ICISTR.GE.2) THEN
            CALL ITODS(IDET,1,-1,LUDIA)
            CALL TODSC_MCLR(DIAG,IDET,-1,LUDIA)
            IDET = 0
          END IF
  900   CONTINUE
  999   CONTINUE
*
 1000 CONTINUE
*
      IF ( ICISTR.GE.2 ) CALL ITODS(-1,1,-1,LUDIA)
*
      RETURN
c Avoid unused argument lines
      IF (.FALSE.) THEN
        CALL Unused_real_array(XA)
        CALL Unused_real_array(SCR)
        CALL Unused_real(PLSIGN)
        CALL Unused_integer(IPRNT)
      END IF
      END
