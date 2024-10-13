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
* Copyright (C) 1992,1994, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE FMAT_CASPT2(FIMO,FAMO,DREF,NBUF,BUF)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "SysDef.fh"
      DIMENSION FIMO(NFIMO),FAMO(NFAMO)
      DIMENSION DREF(NDREF),BUF(NBUF)
      DIMENSION IAD2M(3,36*36)

C COMPUTE THE INACTIVE FOCK MATRIX IN MO BASIS.
C COMPUTE THE ACTIVE FOCK MATRIX IN MO BASIS.
C DEFINITIONS ARE:
C FIMO(P,Q)= H(P,Q)+ SUM( (2*(PQ,KK)-(PK,QK)) OVER INACTIVE K)
C WHERE H(P,Q) IS THE CORE FOCK MATRIX OBTAINED FROM THE
C TRANSFORMATION PART, I.E., THE CI 1-ELECTRON HAMILTONIAN.
C SIMILARLY, FAMO(P,Q)=
C  (1/2)*SUM( (2*(PQ,TU)-(PU,QT))*DREF(T,U) OVER ACTIVE T,U)
C THIS ROUTINE USES DIRECTLY THE TRANSFORMED INTEGRALS TO SET UP
C FIMO AND FAMO.
C CODED 92-12-04 BY MALMQVIST FOR CASPT2, MOLCAS-3 VERSION.


      NDIM2M=(NSYM*(NSYM+1))/2
      IDISK=0
      CALL IDAFILE(LUINTM,2,IAD2M,3*36*36,IDISK)

      IFSTA=1
      DO ISYR=1,NSYM
        NBR=NORB(ISYR)
        IF(NBR.EQ.0) GOTO 120
        NB3=(NBR**2+NBR)/2
        NBNB=NBR**2
        ISYS=ISYR
        IS3RS=(ISYR**2+ISYR)/2
        DO ISYP=1,NSYM
          NIP=NISH(ISYP)
          NOP=NOSH(ISYP)
          NAESP=NAES(ISYP)
          IF(NOP.EQ.0) GOTO 110
          ISYQ=ISYP
          IS3PQ=(ISYP*(ISYP+1))/2
          ISADDR=IS3RS+NDIM2M*(IS3PQ-1)
          IDISK1=IAD2M(1,ISADDR)
          IF(IDISK1.EQ.0) THEN
            WRITE(6,*)' FMAT: COULOMB INTEGRAL BUFFER MISSING!'
            WRITE(6,'(1X,A,4I4)')'SYMMETRY BLOCK:',ISYP,ISYQ,ISYR,ISYS
            CALL ABEND()
          END IF
          DO IT=1,NOP
            DO IU=1,IT
              IF(IT.LE.NIP .AND. IT.EQ.IU) THEN
                CALL DDAFILE(LUINTM,2,BUF,NB3,IDISK1)
C ADD 2* BUFFER INTO CORRECT POSITION OF  FIMO.
                CALL DAXPY_(NB3,2.0D00,BUF,1,FIMO(IFSTA),1)
              ELSE IF(IT.GT.NIP .AND. IT.LE.NOP .AND.
     &                IU.GT.NIP .AND. IU.LE.NOP) THEN
                ITABS=NAESP+(IT-NIP)
                IUABS=NAESP+(IU-NIP)
                ITU=(ITABS*(ITABS-1))/2 + IUABS
                DTU=DREF(ITU)
                IF(IT.EQ.IU) DTU=0.5D0*DTU
                CALL DDAFILE(LUINTM,2,BUF,NB3,IDISK1)
                CALL DAXPY_(NB3,2.0D0*DTU,BUF,1,FAMO(IFSTA),1)
              ELSE
                CALL DDAFILE(LUINTM,0,BUF,NB3,IDISK1)
              END IF
            END DO
          END DO
 110      CONTINUE
        END DO
        IFSTA=IFSTA+NB3
 120    CONTINUE
      END DO

      IFSTA=1
      DO ISYR=1,NSYM
        NBR=NORB(ISYR)
        IF(NBR.EQ.0) GOTO 220
        NB3=(NBR**2+NBR)/2
        NBNB=NBR**2
        ISYS=ISYR
        IS3RS=(ISYR**2+ISYR)/2
        DO ISYP=1,NSYM
          NIP=NISH(ISYP)
          NOP=NOSH(ISYP)
          NAESP=NAES(ISYP)
          IF(NOP.EQ.0) GOTO 210
          ISYQ=ISYP
          IS3PQ=(ISYP*(ISYP+1))/2
          ISADDR=IS3RS+NDIM2M*(IS3PQ-1)
          IDISK2=IAD2M(2,ISADDR)
          IF(IDISK2.EQ.0) THEN
            WRITE(6,*)' FMAT: EXCH-1  INTEGRAL BUFFER MISSING!'
            WRITE(6,'(1X,A,4I4)')'SYMMETRY BLOCK:',ISYP,ISYQ,ISYR,ISYS
            CALL ABEND()
          END IF
          DO IT=1,NOP
            DO IU=1,IT
              IF(IT.LE.NIP .AND. IT.EQ.IU) THEN
                CALL DDAFILE(LUINTM,2,BUF,NBNB,IDISK2)
                CALL TRIANG(NBR,BUF)
C ADD -1* BUFFER INTO CORRECT POSITION OF  FIMO.
                CALL DAXPY_(NB3,-1.0D00,BUF,1,FIMO(IFSTA),1)
              ELSE IF(IT.GT.NIP .AND. IT.LE.NOP .AND.
     &                IU.GT.NIP .AND. IU.LE.NOP) THEN
                ITABS=NAESP+(IT-NIP)
                IUABS=NAESP+(IU-NIP)
                ITU=(ITABS*(ITABS-1))/2 + IUABS
                DTU=DREF(ITU)
                IF(IT.EQ.IU) DTU=0.5D0*DTU
                CALL DDAFILE(LUINTM,2,BUF,NBNB,IDISK2)
                CALL TRIANG(NBR,BUF)
                CALL DAXPY_(NB3,-DTU,BUF,1,FAMO(IFSTA),1)
              ELSE
                CALL DDAFILE(LUINTM,0,BUF,NBNB,IDISK2)
              END IF
            END DO
          END DO
 210      CONTINUE
        END DO
        IFSTA=IFSTA+NB3
 220    CONTINUE
      END DO


      RETURN
      END
