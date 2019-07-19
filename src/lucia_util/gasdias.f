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
*               2015, Lasse Kragh Soerensen                            *
************************************************************************
      SUBROUTINE GASDIAS(    NAEL,   IASTR,    NBEL,   IBSTR,    NORB,
     &                       DIAG,   NSMST,       H,            XB,
     &                            RJ,      RK,   NSSOA,   NSSOB,
     &                      LUDIA,   ECORE,   PSSIGN,   IPRNT,
     &                      NTOOB,  ICISTR,   RJKAA,     I12,   IBLTP,
*
     &                     NBLOCK, IBLKFO,I_AM_OUT,N_ELIMINATED_BATCHES)
*
* Calculate determinant diagonal
* Turbo-ras version
*
* Driven by IBLKFO, May 97
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
* Added the possibility to zero out unwanted diagonal part
* which is needed for highly excited states. Lasse 2015
*
      IMPLICIT REAL*8           (A-H,O-Z)
#include "io_util.fh"
C     REAL * 8  INPROD
*.General input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION H(NORB)
      DIMENSION I_AM_OUT(*)
*. Specific input
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XB(NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION DIAG(*)

      INTEGER IDUM_ARR(1)
*
      NTEST =  0
      NTEST = MAX(NTEST,IPRNT)
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1.0D6
      ELSE
         XADD = 0.0D0
      END IF
*

      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMAT(H,1,NORB,1,NORB)
        WRITE(6,*) ' Core energy ', ECORE
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
        END IF
*
        WRITE(6,*) ' TTSS for Blocks '
        DO IBLOCK = 1, NBLOCK
          WRITE(6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
        END DO
*
        WRITE(6,*) ' I12 = ',I12
      END IF
*
*  Diagonal elements according to Handys formulae
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
      IDET = 0
      ITDET = 0
      IF(LUDIA.NE.0) IDISK(LUDIA)=0
*
      DO IBLK = 1, NBLOCK
* Lasse addition
        I_AM_NOT_WANTED = 0
        DO I = 1, N_ELIMINATED_BATCHES
          IF(I_AM_OUT(I).EQ.IBLK) THEN
            I_AM_NOT_WANTED = 1
            EXIT
          END IF
        END DO
* Lasse addition end
*
        IATP = IBLKFO(1,IBLK)
        IBTP = IBLKFO(2,IBLK)
        IASM = IBLKFO(3,IBLK)
        IBSM = IBLKFO(4,IBLK)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        IDUM_ARR=0
        CALL GETSTR_TOTSM_SPGP(      1,   IATP,   IASM,   NAEL, NASTR1,
     &                           IASTR,   NORB,     0,IDUM_ARR,IDUM_ARR)
        IOFF =  1
        DO IA = 1, NSSOA(IASM,IATP)
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
*. Obtain beta strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP(      2,   IBTP,   IBSM,   NBEL, NBSTR1,
     &                           IBSTR,   NORB,     0,IDUM_ARR,IDUM_ARR)
        IBSTRT = 1
        IBSTOP =  NSSOB(IBSM,IBTP)
        DO IB = IBSTRT,IBSTOP
          IBREL = IB - IBSTRT + 1
*
*. Terms depending only on IB
*
          HB = 0.0D0
          RJBB = 0.0D0
          CALL SETVEC(XB,0.0D0,NORB)
*
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
          IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1
          END IF
          IASTOP = NSSOA(IASM,IATP)
*
          DO IA = IASTRT,IASTOP
            IDET = IDET + 1
            ITDET = ITDET + 1
            X = EB + RJKAA(IA-IOFF+1)
            DO IEL = 1, NAEL
              X = X +XB(IASTR(IEL,IA))
            END DO
* Lasse addition
            IF(I_AM_NOT_WANTED.EQ.0) THEN
              DIAG(IDET) = X
              IF(IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
            ELSE
              DIAG(IDET) = 0.0D0
            END IF
* Lasse addition end
          END DO
*         ^ End of loop over alpha strings|
        END DO
*       ^ End of loop over betastrings
*. Yet a RAS block of the diagonal has been constructed
        IF(ICISTR.GE.2) THEN
          IF(NTEST.GE.100) THEN
            write(6,*) ' number of diagonal elements to disc ',IDET
            CALL WRTMAT(DIAG,1,IDET,1,IDET)
          END IF
          CALL ITODS([IDET],1,-1,LUDIA)
          CALL TODSC(DIAG,IDET,-1,LUDIA)
          IDET = 0
        END IF
      END DO
*        ^ End of loop over blocks

      IF(NTEST.GE.5) WRITE(6,*)
     &' Number of diagonal elements generated (1)',ITDET
*
      IF(NTEST .GE.100 .AND.ICISTR.LE.1 ) THEN
        WRITE(6,*) ' CIDIAGONAL '
        CALL WRTMAT(DIAG(1),1,IDET,1,IDET)
      END IF
*
      IF ( ICISTR.GE.2 ) CALL ITODS([-1],1,-1,LUDIA)
*
      RETURN
      END
