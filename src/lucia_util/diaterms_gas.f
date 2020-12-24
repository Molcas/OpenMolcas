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
* Copyright (C) 1995, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE DIATERMS_GAS(   NAEL,  IASTR,   NBEL,  IBSTR,   NORB,
     &                            VEC,  NSMST,      H,    IDC,
     &                             XB,             RJ,     RK,  NSSOA,
     &                          NSSOB,  ECORE,   LUIN,  LUOUT,  IPRNT,
     &                          NTOOB,  RJKAA,    I12, IBLOCK, NBLOCK,
*
     &                          ITASK, FACTOR,  I0CHK,  I0BLK)
*
* Terms from diagonal to specific blocks
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
*
* Calculate determinant diagonal
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, July 1995, GAS version
*
* I12 = 1 => only one-body part
*     = 2 =>      one+two-body part
*
      IMPLICIT REAL*8           (A-H,O-Z)
#include "io_util.fh"
*.General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION H(NORB)
      DIMENSION IBLOCK(8,*)
*.
      INTEGER I0BLK(*)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XB(NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION VEC (*)

      INTEGER IDUM_ARR(1)
*
      NTEST =  00
      NTEST = MAX(NTEST,IPRNT)
C?    WRITE(6,*) ' NTEST = ',NTEST
*
      IF(LUIN.GT.0) IDISK(LUIN)=0
      IF(LUOUT.GT.0) IDISK(LUOUT)=0

      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' DIATERMS_GAS in action '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' LUIN,LUOUT = ', LUIN,LUOUT
        WRITE(6,*) ' NBLOCK =', NBLOCK
        WRITE(6,*) ' I0CHK = ', I0CHK
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
          WRITE(6,*) ' I12 and ITASK = ', I12,ITASK
        END IF
      WRITE(6,*) ' FACTOR = ',FACTOR
      END IF
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
      IF(I12.EQ.2)
     &CALL VECSUM(RK,RK,RJ,-1.0D0,+1.0D0,NTOOB **2)
*
      ITDET = 0
      IDET = 0
      DO JBLOCK = 1, NBLOCK
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
        IOFF = IBLOCK(6,JBLOCK)
        IF(NTEST.GE.20) THEN
         WRITE(6,*) ' Block in action : IATP IBTP IASM IBSM ',
     &               IATP,IBTP,IASM,IBSM
        END IF
*
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
*
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM_ARR = 0
        CALL GETSTR_TOTSM_SPGP(      1,   IATP,   IASM,   NAEL, NASTR1,
     &                           IASTR,   NORB,     0,IDUM_ARR,IDUM_ARR)
        IF(NTEST.GE.1000) THEN
          write(6,*) ' After GETSTR for A strings '
          WRITE(6,*) ' alpha strings obtained '
          NAST = NSSOA(IASM,IATP)
          CALL IWRTMA(IASTR,NAEL,NAST,NAEL,NAST)
        END IF
*
        IOFF =  1
        NIA = NSSOA(IASM,IATP)
        DO IA = 1 ,NSSOA(IASM,IATP)
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
          RJKAA(IA-IOFF+1) = EAA
        END DO
*. Obtain alpha strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP(      2,   IBTP,   IBSM,   NBEL, NBSTR1,
     &                           IBSTR,   NORB,     0,IDUM_ARR,IDUM_ARR)
        NIB =  NSSOB(IBSM,IBTP)
*
        IMZERO=0
        IF(LUIN.GT.0) THEN
          CALL IDAFILE(LUIN,2,IDUM_ARR,1,IDISK(LUIN))
          LDET=IDUM_ARR(1)
          CALL IDAFILE(LUIN,2,IDUM_ARR,1,IDISK(LUIN))
          IDET = 0
          CALL FRMDSC(   VEC(1),     LDET,       -1,     LUIN,   IMZERO,
     &                  IAMPACK)
        END IF
*
        IF(I0CHK.EQ.1) THEN
          IMZERO = I0BLK(JBLOCK)
          IF(IMZERO.EQ.1) THEN
*.Update offset to next block
            IF(IPACK.EQ.1.AND.IATP.EQ.IBTP) THEN
              IDET = IDET + NIA*(NIA+1)/2
            ELSE
              IDET = IDET + NIA*NIB
            END IF
          END IF
        END IF
C?      WRITE(6,*) ' DIATERMS_GAS : I0CHK,JBLOCK IMZERO',
C?   &  I0CHK,JBLOCK,IMZERO
*
        IF(IMZERO.NE.1) THEN
*. Calculate ...
*
        DO IB = 1 ,NIB
*
*. Terms depending only on IB
*
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
*
          IF(IPACK.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1
          END IF
*
          IASTOP = NSSOA(IASM,IATP)
          DO IA = IASTRT,IASTOP
            IDET = IDET + 1
            ITDET = ITDET + 1
            X = EB + RJKAA(IA-IOFF+1)
            DO IEL = 1, NAEL
              X = X +XB(IASTR(IEL,IA))
            END DO
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
            IF(ITASK.EQ.1) THEN
              IF(ABS(X+FACTOR) .GT. 1.0D-10) THEN
                VEC(IDET) = VEC(IDET)/(X+FACTOR)
              ELSE
                VEC(IDET) = 0.0D0
              END IF
            ELSE
              VEC(IDET) = VEC(IDET)*(X+FACTOR)
            END IF
C?         write(6,*) ' IDET,X,VEC(IDET) ', IDET,X,VEC(IDET)
          END DO
        END DO
        END IF
*
        IF(LUOUT.GT.0) THEN
          CALL ITODS([LDET],1,-1,LUOUT)
          CALL TODSC(VEC,LDET,-1,LUOUT)
C?        WRITE(6,*) ' Number of elements transferred to DISC ',
C?   &    LDET
          IDET = 0
        END IF
*
      END IF
      END DO
*
      IF(LUOUT.GT.0) THEN
       CALL ITODS([-1],1,-1,LUOUT)
      END IF
*
C?    WRITE(6,*) ' Mission DIATERMS finished '
*
      RETURN
      END
