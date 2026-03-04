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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE TRDACT(IVEC,JVEC,DTU)
      use gugx, only: SGS
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_global, only: LUCIEX, IDTCEX
      use caspt2_module, only: nAshT, iSCF, jState, nConf, nSym, STSym,
     &                         iASym, nAes, nAshT, nAsh
      use gugx, only: MxLev
      use pt2_guga, only: MxCI
      use constants, only: Zero, One, Two
      use definitions, only: iwp, wp
      IMPLICIT NONE

      INTEGER(kind=iwp) IVEC, JVEC
      REAL(kind=wp) DTU(NASHT,NASHT)
C Local array:
      INTEGER(kind=iwp) IATOG(MXLEV)

      REAL(kind=wp), ALLOCATABLE :: TRDOP1(:), TRDOP2(:), TRDOP3(:)
      INTEGER(kind=iwp) :: NOP1, NOP2, NOP3

      REAL(kind=wp), ALLOCATABLE :: TRDTMP(:), TRDCI(:), TRDSGM(:)

      INTEGER(kind=iwp) :: I, ID
      INTEGER(kind=iwp) :: ISYM, ISYMT
      INTEGER(kind=iwp) :: ITABS, ITLEV, IU, IUABS, IULEV
      REAL(kind=wp) :: OP0, OCCNUM, SCP, DDOT_
      Integer(kind=iwp) :: nLev

      nLev = SGS%nLev

C Add to the active-active block of transition density matrix,
C    D(t,u) = Add <IVEC| E(t,u) |JVEC> = <0| W1T E(t,u) W2 |0>
C where t,u are active indices. IVEC and JVEC are integer labels,
C denoting sets of coefficients stored on LUSOLV. These are assumed
C to be contravariant representations of the wave operators W1 and W2,
C in the notation of the comments.


C (1): Compute a representation of the operator PCAS*W1T*W2
      NOP1=NASHT**2
      NOP2=(NOP1*(NOP1+1))/2
      NOP3=(NOP2*(NOP1+2))/3
      CALL MMA_ALLOCATE(TRDOP1,NOP1)
      CALL MMA_ALLOCATE(TRDOP2,NOP2)
      CALL MMA_ALLOCATE(TRDOP3,NOP3)
      CALL MKWWOP(IVEC,JVEC,OP0,TRDOP1,NOP2,TRDOP2,NOP3,TRDOP3)


C (2): Compute the state vector |Temp> = (PCAS*W1T*W2) |0>
C First modify the coefficients, see subroutine MODOP.
      CALL MODOP(TRDOP1,NOP2,TRDOP2,NOP3,TRDOP3)
      CALL MMA_ALLOCATE(TRDTMP,NCONF)
      CALL MMA_ALLOCATE(TRDCI,NCONF)
      IF(ISCF.EQ.0) THEN
        ID=IDTCEX(JSTATE)
        CALL DDAFILE(LUCIEX,2,TRDCI,NCONF,ID)
      ELSE
        TRDCI(1)=One
      END IF
      TRDTMP(:)=Zero
      CALL HAM3(OP0,TRDOP1,NOP2,TRDOP2,NOP3,TRDOP3,STSYM,
     &          TRDCI,TRDTMP,nCONF)
C No more need for the operators:
      CALL MMA_DEALLOCATE(TRDOP1)
      CALL MMA_DEALLOCATE(TRDOP2)
      CALL MMA_DEALLOCATE(TRDOP3)


C (3) compute <0| E(t,u) W1T W2 |0> as <Sigma_ut| Temp>, and add to DTU

      IF(ISCF.EQ.0) THEN
C Create reorder table giving the GUGA level, i.e. CI-coupling
C ordinal number of each active orbital.
        ITABS=0
        DO ISYM=1,NSYM
          DO I=1,NLEV
            IF(SGS%ISM(I).EQ.ISYM) THEN
              ITABS=ITABS+1
              IATOG(ITABS)=I
            END IF
          END DO
        END DO
        CALL MMA_ALLOCATE(TRDSGM,MXCI)
        DO ITABS=1,NASHT
          ISYMT=IASYM(ITABS)
          ITLEV=IATOG(ITABS)
          DO IU=1,NASH(ISYMT)
            IUABS=NAES(ISYMT)+IU
            IULEV=IATOG(IUABS)
            CALL GETSGM2(IULEV,ITLEV,STSYM,TRDCI,NCONF,TRDSGM,NCONF)
            SCP=DDOT_(NCONF,TRDSGM,1,TRDTMP,1)
            DTU(ITABS,IUABS)=DTU(ITABS,IUABS)+SCP
          END DO
        END DO
        CALL MMA_DEALLOCATE(TRDSGM)
      ELSE
        OCCNUM=Two
        IF(ISCF.EQ.2) OCCNUM=One
        DO ITABS=1,NASHT
          DTU(ITABS,ITABS)=DTU(ITABS,ITABS)+OCCNUM
        END DO
      END IF
C No more need for CI array
      CALL MMA_DEALLOCATE(TRDCI)
C No more need for the TMP state vector
      CALL MMA_DEALLOCATE(TRDTMP)

C (4): Add the correction <0| [W1T,E(tu)] W2 |0>.
      CALL COMMWEW(IVEC,JVEC,DTU)

      END SUBROUTINE TRDACT
