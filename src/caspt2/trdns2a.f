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
      SUBROUTINE TRDNS2A(IVEC,JVEC,DPT2,NDPT2)

      use definitions, only: iwp, wp
      use constants, only: Zero, Two
      use caspt2_global, only:iPrGlb
      use caspt2_global, only: DREF
      use PrintLevel, only: verbose
      use caspt2_module, only: nActEl, nAshT, nSym, nInDep, nISup,
     &                         nIsh, nAsh, nOrb, nAES
      IMPLICIT None


      integer(kind=iwp), intent(in):: IVEC, JVEC, NDPT2
      real(kind=wp), intent(inout):: DPT2(NDPT2)

      integer(kind=iwp) ::
     &               NACTD(13)=[1, 2, 2,-1, 0, 1, 1,-2,-2,-1,-1, 0, 0]
      real(kind=wp) COEF1, COEF2, D, DR, OVL
      integer(kind=iwp) ICASE, IOFDPT, ISYM, IT, ITABS, ITQ, ITU, IU,
     &                  IUABS, IUQ, IUT, lVec1, lVec2, NA, NADIFF,
     &                  NAHOLE, NI, NIN, NIS, NO, nVec
      real(kind=wp), External:: RHS_DDOT

C Add to the diagonal blocks of transition density matrix,
C    DPT2(p,q) = Add <IVEC| E(p,q) |JVEC>,
C where p,q are active indices. Compare TRDNS2D.
C The present solution gives just a reasonable approximation,
C with correct trace.
      IF ( IPRGLB.GE.VERBOSE ) THEN
      Call WarningMessage(1,'Computing approximated density.')
      WRITE(6,*)' The active/active submatrices of the density'
      WRITE(6,*)' matrix is roughly approximated only.'
      END IF

      COEF1=Zero
      COEF2=Zero
      NAHOLE=2*NASHT-NACTEL
      DO ICASE=1,13
        NADIFF=NACTD(ICASE)
        IF(NACTEL+NADIFF.LT.0) Cycle
        IF(NAHOLE-NADIFF.LT.0) Cycle
        OVL=Zero
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) Cycle
          NIS=NISUP(ISYM,ICASE)
          NVEC=NIN*NIS
          IF(NVEC.EQ.0) Cycle
          CALL RHS_ALLO(NIN,NIS,LVEC1)
          CALL RHS_ALLO(NIN,NIS,LVEC2)
          CALL RHS_READ_SR (LVEC1,iCASE,iSYM,IVEC)
          CALL RHS_READ_SR (LVEC2,iCASE,iSYM,JVEC)
          OVL=OVL+RHS_DDOT(NIN,NIS,LVEC1,LVEC2)
          CALL RHS_FREE(LVEC1)
          CALL RHS_FREE(LVEC2)
        End Do
        IF(NADIFF.GT.0) THEN
          COEF1=COEF1+OVL*DBLE(NADIFF)/DBLE(MAX(1,NAHOLE))
          COEF2=COEF2+OVL*DBLE(NAHOLE-NADIFF)/DBLE(MAX(1,NAHOLE))
        ELSE
          COEF2=COEF2+OVL*DBLE(NACTEL+NADIFF)/DBLE(MAX(1,NACTEL))
        END IF
      End Do

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
            DR=DREF((ITABS*(ITABS-1))/2+IUABS)
            D=COEF2*DR
            IF(IT.EQ.IU) D=D+Two*COEF1
            ITU=ITQ+NO*(IUQ-1)
            IUT=IUQ+NO*(ITQ-1)
            DPT2(IOFDPT+ITU)=DPT2(IOFDPT+ITU)+D
            DPT2(IOFDPT+IUT)=DPT2(IOFDPT+ITU)
          END DO
        END DO
        IOFDPT=IOFDPT+NO**2
      END DO

      END SUBROUTINE TRDNS2A
