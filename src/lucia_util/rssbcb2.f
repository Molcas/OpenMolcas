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
* Copyright (C) 1991, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE RSSBCB2(    IASM,    IATP,      IBSM,    IBTP,
     &                        JASM,    JATP,    JBSM,    JBTP,
     &                       NGAS,    IAOC,    IBOC,    JAOC,    JBOC,
     &                       NAEL,    NBEL,  IJAGRP,  IJBGRP,      SB,
     &                         CB,   JDOH2,   ADSXA,    STSTSX,
*
     &                     DXSTST,  STSTDX,  SXDXSX,  NOBPTS,  IOBPTS,
     &                    MXPNGAS,   ITSOB,    MAXI,    MAXK,    SSCR,
     &                       CSCR,      I1,    XI1S,      I2,    XI2S,
     &                       XINT,      C2,   NSMOB,   NSMST,   NSMSX,
     &                      NSMDX,     NIA,     NIB,     NJA,     NJB,
*
     &                     MXPOBS,     IDC,         CJRES,
     &                      SIRES,      I3,    XI3S,      I4,    XI4S,
     &                     MXSXBL,  MXSXST,   MOCAA,
     &                      IPRNT,   IHAPR,
     &                       SCLFAC, IUSE_PH,  IPHGAS,
*
     &                   I_RES_AB,
     &                     XINT2)
*
* SUBROUTINE RSSBCB2 --> 82
*
*
* Contributions to sigma block (iasm iatp, ibsm ibtp ) from
* C block (jasm jatp , jbsm, jbtp)
*
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha strings in sigma
* IBSM,IBTP : Symmetry and type of beta  strings in sigma
* JASM,JATP : Symmetry and type of alpha strings in C
* JBSM,JBTP : Symmetry and type of beta  strings in C
* NGAS      : Number of active spaces in calculation
* IAOC,IBOC : Number of electrons in each AS for sigma supergroups
* JAOC,JBOC : Number of electrons in each AS for C     supergroups
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* IDOH2 : = 0 => no two electron operator
* IDOH2 : = 1 =>    two electron operator
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
*
* IHAPR : .ne. 0 implies thatt the exact Hamiltonian shoulf not be uses
* In the case IPTSPC and JPTSPC defined the perturbation spaces
* a nonvanishing perturbation is allowed inside each subspace.
* The actual type of approximate Hamiltonian in each subspace is defined by
* IHFORM
* NNSEL2E : Only selected 2e terms will be included
* ISEL2E : orbital spaces in which 2e terms are included
*          (Currently : all indeces identical )
*
* ======
* Output
* ======
* SB : fresh sigma block
*
* =======
* Scratch
* =======
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* C2 : Must hold largest STT block of sigma or C
*
* XINT : Scratch space for integrals.
*
* Jeppe Olsen , Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "timers.fh"
      INTEGER  ADSXA(*),STSTSX(*),DXSTST(*),STSTDX(*),SXDXSX(*)
*. Output
      DIMENSION CB(*),SB(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION C2(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION IBLOCK(8)
      DIMENSION IPHGAS(*)
      DIMENSION XINT(*),XINT2(*)

      DIMENSION ITSOB(*),NOBPTS(*),IOBPTS(*)
      DIMENSION IAOC(*),JAOC(*),IBOC(*),JBOC(*)
*. For H(apr)
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' ==============================='
        WRITE(6,*) ' RSSBCB2 :  C block (transposed)'
        WRITE(6,*) ' ================================'
        CALL WRTMAT(CB,NJB,NJA,NJB,NJA)
        WRITE(6,*) ' ======================================='
        WRITE(6,*) ' RSSBCB2 : Initial  S block(transposed) '
        WRITE(6,*) ' ======================================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
        WRITE(6,*) ' Overall scalefactor ',SCLFAC
        WRITE(6,*) ' IHAPR,JDOH2 = ', IHAPR,JDOH2
        WRITE(6,*) ' IUSE_PH,I_RES_AB = ', IUSE_PH,I_RES_AB
      END IF
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' IAOC and IBOC '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' JAOC and JBOC  : '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' IASM IATP JASM JATP ', IASM,IATP,JASM,JATP
        WRITE(6,*) ' IBSM IBTP JBSM JBTP ', IBSM,IBTP,JBSM,JBTP
        WRITE(6,*) ' NAEL NBEL ', NAEL, NBEL
      END IF
* Should the corresponding Hamiltonian matrix block be
* calculated exactly or approximately
c      IF(IHAPR.NE.0) THEN
c        CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
c     &               IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
c     &               IMZERO,IDIAG)
c        IF(NTEST.GE. 20)
c     &  WRITE(6,*) ' RSSBCBN : ', NNSEL2E, ISEL2E(1)
c        NSEL2E = NNSEL2E
c        IF(IMZERO.NE.0) GOTO 9999
c      ELSE
*. Operator specified by input
        IAPRLEV =-1
        IDOH2 = JDOH2
        IDIAG = 0
        NSEL2E = 0
c      END IF
      IF(NTEST.GE. 20)
     &WRITE(6,*) ' IHAPR, IDIAG IDOH2 ' , IHAPR,IDIAG, IDOH2
*
*
      IF(IDC.EQ.2.AND.IATP.EQ.IBTP.AND.IASM.EQ.IBSM .AND.
     &   I_RES_AB.EQ.0.AND.JASM.EQ.JBSM.AND.JATP.EQ.JBTP) THEN
*. Diagonal sigma block, use alpha-beta symmetry to reduce computations.
        IUSEAB = 1
      ELSE
        IUSEAB = 0
      END IF
*
      IF(IDIAG.EQ.0) THEN
*
* Calculate block exactly
*
      IF(I_RES_AB.NE.1.AND.IUSEAB.EQ.0.
     &   AND.IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
* Sigma beta beta contribution
* =============================
*
* Sigma aa(IA,IB) = sum(i.gt.k,j.gt.l)<IB!Eb(ij)Eb(kl)!JB>
*                 * ((ij!kl)-(il!kj)) C(IA,JB)
*                 + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)
*.One electron part
          CALL TRPMT3(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)

        IF(NBEL.GE.0) THEN
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB before RSBB1E'
             call wrtmat(sb,nia,nib,nia,nib)
           END IF
          IF(NTEST.GE.101)
     &    WRITE(6,*) ' I am going to call RSBB1E'
          CALL TIMING(CPU0,CPU,WALL0,WALL)
          CALL RSBB1E_LUCIA(   IBSM,   IBTP,   JBSM,   JBTP, IJBGRP,
     &                          NIA,   NGAS,   IBOC,   JBOC,     SB,
     &                           CB,  ADSXA, STSTSX, NOBPTS,   MAXI,
     &                         MAXK,   SSCR,   CSCR,     I1,   XI1S,
     &                           I2,   XI2S,   XINT,  NSMOB,  NSMST,
*
     &                        NSMSX,  MOCAA, MXSXST,  MOCAA, SCLFAC,
     &                      IUSE_PH, IPHGAS,  NTEST)
          CALL TIMING(CPU1,CPU,WALL1,WALL)
          TSIGMA(1)=TSIGMA(1)+(WALL1-WALL0)
*
* CALL RSBB1E_LUCIA --> 33
*
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB after RSBB1E'
             call wrtmat(sb,nib,nia,nib,nia)
           END IF
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' first element of SB after RSBB1E',
     &       SB(1)
           END IF

        END IF
        IF(IDOH2.NE.0.AND.NBEL.GE.0) THEN
*. Two electron part
          IF(NTEST.GE.101)
     &    WRITE(6,*) ' I am going to call RSBB2A'
          CALL TIMING(CPU0,CPU,WALL0,WALL)
          CALL RSBB2A_LUCIA(   IBSM,   IBTP,   JBSM,   JBTP, IJBGRP,
     &                          NIA,   NIB,    NGAS,   IBOC,   JBOC,
     &                           SB,     CB,  ADSXA,  STSTDX,
     &                       SXDXSX, NOBPTS,
     &                         MAXI,   MAXK,   SSCR,   CSCR,     I1,
*
     &                         XI1S,    XINT,  NSMOB,
     &                        NSMST,    NSMDX,
     &                         MXSXST, MXSXBL,  MOCAA, SCLFAC,
     &                        NTEST,   IPHGAS)
          CALL TIMING(CPU1,CPU,WALL1,WALL)
          TSIGMA(2)=TSIGMA(2)+(WALL1-WALL0)
*
*
* CALL RSBB2A_LUCIA --> 46
*
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB after RSBB2a'
             call wrtmat(sb,nib,nia,nib,nia)
           END IF
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' first element of SB after RSBB1E',
     &       SB(1)
           END IF
        END IF
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
      END IF
*
* =============================
* Sigma alpha beta contribution
* =============================
*
      IF(IDOH2.NE.0.AND.NAEL.GE.0.AND.NBEL.GE.0) THEN
        IF(NTEST.GE.101)
     &  WRITE(6,*) ' I am going to call RSBB2B'
        IIITRNS = 1
        IF(IIITRNS.EQ.1) THEN
*. Call advice routine
C     ADVICE_SIGMA(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,LADVICE)
           CALL ADVICE_SIGMA(   IAOC,   IBOC,   JAOC,   JBOC, LADVICE)
*. LADVICE = 2 => implies transpose
           IF(LADVICE.EQ.2) THEN
             JJJTRNS = 1
           ELSE
             JJJTRNS = 0
           END IF
        END IF
*
C       WRITE(6,*) ' IUSE_PA = ', IUSE_PA
*
        IF (JJJTRNS.EQ.0) THEN
c          IF( IUSE_PA.EQ.0 ) THEN
          CALL TIMING(CPU0,CPU,WALL0,WALL)
          CALL RSBB2BN_LUCIA(   IASM,   IATP,   IBSM,   IBTP,    NIA,
     &                           NIB,   JASM,   JATP,   JBSM,   JBTP,
     &                           NJA,    NJB, IJAGRP, IJBGRP,   NGAS,
     &                          IAOC,   IBOC,   JAOC,   JBOC,     SB,
     &                            CB,  ADSXA, STSTSX,MXPNGAS, NOBPTS,
*
     &                          MAXK,   SSCR,   CSCR,     I1,   XI1S,
     &                            I2,   XI2S,     I3,   XI3S,     I4,
     &                          XI4S,   XINT,  NSMOB,  NSMST,  NSMSX,
     &                         NSMDX, MXPOBS, IUSEAB,  CJRES,  SIRES,
     &                        SCLFAC,  NTEST,      0,    [0],IUSE_PH,
*
     &                        IPHGAS,  XINT2)
          CALL TIMING(CPU1,CPU,WALL1,WALL)
          TSIGMA(3)=TSIGMA(3)+(WALL1-WALL0)
*
* CALL RSBB2BN_LUCIA --> 52
*
c          ELSE IF ( IUSE_PA.EQ.1 ) THEN
cC         WRITE(6,*) ' RSBB2BVN will be called '
c          CALL RSBB2BVN(IASM,IATP,IBSM,IBTP,NIA,NIB,
c     &         JASM,JATP,JBSM,JBTP,NJA,NJB,
c     &         IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
c     &         SB,CB,ADSXA,STSTSX,MXPNGAS,
c     &         NOBPTS,MAXK,
c     &         SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c     &         XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
c     &         IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
c     &         CJPA,SIPA)
c          END IF
*
         ELSE IF ( JJJTRNS.EQ.1) THEN
*. well lets give the transpose routine some more practice : Transpose back
          CALL TRPMT3(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
*
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
C         WRITE(6,*) ' RSSBCB2 : Transpose path choosen'
*
c          IF(IUSE_PA.EQ.0) THEN
*. No division into active/passive
          CALL TIMING(CPU0,CPU,WALL0,WALL)
          CALL RSBB2BN_LUCIA(   IBSM,   IBTP,   IASM,   IATP,    NIB,
     &                           NIA,   JBSM,   JBTP,   JASM,   JATP,
     &                           NJB,    NJA, IJBGRP, IJAGRP,   NGAS,
     &                          IBOC,   IAOC,   JBOC,   JAOC,     SB,
     &                            CB,  ADSXA, STSTSX,MXPNGAS, NOBPTS,
*
     &                          MAXK,   SSCR,   CSCR,     I1,   XI1S,
     &                            I2,   XI2S,     I3,   XI3S,     I4,
     &                          XI4S,   XINT,  NSMOB,  NSMST,  NSMSX,
     &                         NSMDX, MXPOBS, IUSEAB,  CJRES,  SIRES,
     &                        SCLFAC,  NTEST,      0,    [0],IUSE_PH,
*
     &                        IPHGAS,  XINT2)
          CALL TIMING(CPU1,CPU,WALL1,WALL)
          TSIGMA(3)=TSIGMA(3)+(WALL1-WALL0)
*
* CALL RSBB2BN_LUCIA --> 52
*
c          ELSE
c*. Divide into active/passive
c          CALL RSBB2BVN(IBSM,IBTP,IASM,IATP,NIB,NIA,
c     &         JBSM,JBTP,JASM,JATP,NJB,NJA,
c     &         IJBGRP,IJAGRP,NGAS,IBOC,IAOC,JBOC,JAOC,
c     &         SB,CB,ADSXA,STSTSX,MXPNGAS,
c     &         NOBPTS,MAXK,
c     &         SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c     &         XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
c     &         IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
c     &         CJPA,SIPA)
c           END IF

*. Transpose ( To compensate later transposition )
          CALL TRPMT3(SB,NIA,NIB,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJA,NJB,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
        END IF
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB after RSBB2B, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB after RSBB2b'
             call wrtmat(sb,nia,nib,nia,nib)
           END IF
      END IF
*
* =============================
* Sigma alpha alpha contribution
* =============================
*
      IF(I_RES_AB.NE.-1.AND.
     &   NAEL.GE.0.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* alpha single excitation
*
        IF(NTEST.GE.101)
     &  WRITE(6,*) ' I am going to call RSBB1E (last time )'
        CALL TIMING(CPU0,CPU,WALL0,WALL)
        CALL RSBB1E_LUCIA(    IASM,    IATP,    JASM,    JATP,  IJAGRP,
     &                         NIB,    NGAS,    IAOC,    JAOC,      SB,
     &                          CB,   ADSXA,  STSTSX,  NOBPTS,    MAXI,
     &                        MAXK,    SSCR,    CSCR,      I1,    XI1S,
     &                          I2,    XI2S,    XINT,   NSMOB,   NSMST,
*
     &                       NSMSX,   MOCAA,  MXSXST,   MOCAA,  SCLFAC,
     &                     IUSE_PH,  IPHGAS,   NTEST)
        CALL TIMING(CPU1,CPU,WALL1,WALL)
        TSIGMA(1)=TSIGMA(1)+(WALL1-WALL0)
*
* CALL RSBB1E_LUCIA --> 33
*
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB transposed after RSBB1, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB transposed  after RSBB1E'
             call wrtmat(SB,nib,nia,nib,nia)
           END IF
*
* alpha double excitation
*
        IF(IDOH2.NE.0.AND.NAEL.GE.0) THEN
          IF(NTEST.GE.101)
     &    WRITE(6,*) ' I am going to call RSBB2A (last time )'
          CALL TIMING(CPU0,CPU,WALL0,WALL)
          CALL RSBB2A_LUCIA(   IASM,   IATP,   JASM,   JATP, IJAGRP,
     &                          NIB,   NIA,    NGAS,   IAOC,   JAOC,
     &                           SB,     CB,  ADSXA,  STSTDX,
     &                       SXDXSX, NOBPTS,
     &                         MAXI,   MAXK,   SSCR,   CSCR,     I1,
*
     &                         XI1S,       XINT,  NSMOB,
     &                        NSMST,   NSMDX,
     &                         MXSXST, MXSXBL,  MOCAA, SCLFAC,
     &                        NTEST,   IPHGAS)
          CALL TIMING(CPU1,CPU,WALL1,WALL)
          TSIGMA(2)=TSIGMA(2)+(WALL1-WALL0)
*
*
* CALL RSBB2A_LUCIA --> 46
*
        END IF
*
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB transposed after RSBB2A, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
        IF(NTEST.GE.500) THEN
          WRITE(6,*) ' SB after RSBB2A'
          call wrtmat(sb,nia,nib,nia,nib)
        END IF
      END IF
*
      ELSE IF (IDIAG.EQ.1) THEN
*
*. Diagonal approxiation (IDIAG = 1)
*  or complete orbital conserving part of Ham (IH_OCC_CONS = 1)
*
       IBLOCK(1) = IATP
       IBLOCK(2) = IBTP
       IBLOCK(3) = IASM
       IBLOCK(4) = IBSM
       IBLOCK(5) = 1
       IBLOCK(6) = 1
       IF(IDOH2.EQ.0) THEN
         I12 = 1
       ELSE
         I12 = 2
       END IF
C?     WRITE(6,*) ' IDOH2, I12 ', IDOH2,I12
       ITASK = 2
       FACTOR = 0.0D0
*. Well, we are not using transposed matrices here so
       CALL TRPMT3(CB,NJB,NJA,C2)
       CALL COPVEC(C2,CB,NJA*NJB)

       IF(IATP.EQ.JATP.AND.IBTP.EQ.JBTP.AND.
     &    IASM.EQ.JASM.AND.IBSM.EQ.JBSM) THEN
C?       WRITE(6,*) ' DIATERM2_GAS will be called '
         CALL COPVEC(CB,C2,NJA*NJB)
*. Input is in det basis
         IIDC = 1
         CALL DIATERM2_GAS(  FACTOR,   ITASK,      C2,       1,  IBLOCK,
     &                            1,       0,     I12,    IIDC)
       ELSE
         ZERO = 0.0D0
         CALL SETVEC(C2,ZERO,NIA*NIB)
       END IF
*. Remaining occupation conserving operator
c       IF(IH_OCC_CONS.EQ.1) THEN
c         CALL HCONFDIA_BBM(NAEL,NBEL,IJAGRP,IJBGRP,
c     &        IASM,IATP,IAOC,NIA,IBSM,IBTP,IBOC,NIB,
c     &        JASM,JATP,JAOC,NJA,JBSM,JBTP,JBOC,NJB,XINT,CB,C2)
c       END IF

       IF(IUSEAB.EQ.0) THEN
         FACTOR = 1.0D0*SCLFAC
       ELSE
         FACTOR = 0.5D0*SCLFAC
       END IF
C           MAT_P_MATT(A,B,NR,NC,COEF)
       CALL MAT_P_MATT(SB,C2,NIB,NIA,FACTOR)
       CALL TRPMT3(CB,NJA,NJB,C2)
       CALL COPVEC(C2,CB,NJA*NJB)
      END IF
*
c 9999 CONTINUE
c      IF(IHAPR.NE.0) THEN
c*. Clean up
c        CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
c     &               IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
c     &               IMZERO,IDIAG)
c      END IF
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' ==================================='
        WRITE(6,*) ' RSSBCB : Final S block (transposed)'
        WRITE(6,*) ' ==================================='
        CALL WRTMAT(SB,NIB,NIA,NIB,NIA)
      END IF
      NTESTO = NTEST
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(DXSTST)
        CALL Unused_integer_array(IOBPTS)
        CALL Unused_integer_array(ITSOB)
      END IF
      END
