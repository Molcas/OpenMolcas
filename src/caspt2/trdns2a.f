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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*

      SUBROUTINE TRDNS2A(IVEC,JVEC,DPT2)

      IMPLICIT REAL*8 (A-H,O-Z)


#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "sigma.fh"

      DIMENSION DPT2(*)
      DIMENSION NACTD(13)
      DATA NACTD / 1, 2, 2,-1, 0, 1, 1,-2,-2,-1,-1, 0, 0 /

C Add to the diagonal blocks of transition density matrix,
C    DPT2(p,q) = Add <IVEC| E(p,q) |JVEC>,
C where p,q are active indices. Compare TRDNS2D.
C The present solution gives just a reasonable approximation,
C with correct trace.
      CALL QENTER('TRDNS2A')
      IF ( IPRGLB.GE.VERBOSE ) THEN
      Call WarningMessage(1,'Computing approximated density.')
      WRITE(6,*)' The active/active submatrices of the density'
      WRITE(6,*)' matrix is roughly approximated only.'
      END IF

      COEF1=0.0D0
      COEF2=0.0D0
      NAHOLE=2*NASHT-NACTEL
      DO 101 ICASE=1,13
        NADIFF=NACTD(ICASE)
        IF(NACTEL+NADIFF.LT.0) GOTO 101
        IF(NAHOLE-NADIFF.LT.0) GOTO 101
        OVL=0.0D0
        DO 100 ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 100
          NIS=NISUP(ISYM,ICASE)
          NVEC=NIN*NIS
          IF(NVEC.EQ.0) GOTO 100
          !CALL GETMEM('VEC1','ALLO','REAL',LVEC1,NVEC)
          !CALL GETMEM('VEC2','ALLO','REAL',LVEC2,NVEC)
          CALL RHS_ALLO(NIN,NIS,LVEC1)
          CALL RHS_ALLO(NIN,NIS,LVEC2)
          CALL RHS_READ_SR (LVEC1,iCASE,iSYM,IVEC)
          CALL RHS_READ_SR (LVEC2,iCASE,iSYM,JVEC)
          !OVL=OVL+DDOT_(NVEC,WORK(LVEC1),1,WORK(LVEC2),1)
          OVL=OVL+RHS_DDOT(NIN,NIS,LVEC1,LVEC2)
          !CALL GETMEM('VEC1','FREE','REAL',LVEC1,NVEC)
          !CALL GETMEM('VEC2','FREE','REAL',LVEC2,NVEC)
          CALL RHS_FREE(NIN,NIS,LVEC1)
          CALL RHS_FREE(NIN,NIS,LVEC2)
 100    CONTINUE
        IF(NADIFF.GT.0) THEN
          COEF1=COEF1+OVL*DBLE(NADIFF)/DBLE(MAX(1,NAHOLE))
          COEF2=COEF2+OVL*DBLE(NAHOLE-NADIFF)/DBLE(MAX(1,NAHOLE))
        ELSE
          COEF2=COEF2+OVL*DBLE(NACTEL+NADIFF)/DBLE(MAX(1,NACTEL))
        END IF
 101  CONTINUE

      IOFDPT=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NO=NORB(ISYM)
        DO IT=1,NA
          ITQ=NI+IT
          ITABS=NAES(ISYM)+IT
          DO IU=1,IT
            IUQ=NI+IU
            IUABS=NAES(ISYM)+IU
            DR=WORK(LDREF-1+(ITABS*(ITABS-1))/2+IUABS)
            D=COEF2*DR
            IF(IT.EQ.IU) D=D+2.0D0*COEF1
            ITU=ITQ+NO*(IUQ-1)
            IUT=IUQ+NO*(ITQ-1)
            DPT2(IOFDPT+ITU)=DPT2(IOFDPT+ITU)+D
            DPT2(IOFDPT+IUT)=DPT2(IOFDPT+ITU)
          END DO
        END DO
        IOFDPT=IOFDPT+NO**2
      END DO

      CALL QEXIT('TRDNS2A')
      RETURN
      END
