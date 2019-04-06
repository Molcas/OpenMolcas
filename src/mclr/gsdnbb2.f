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
      SUBROUTINE GSDNBB2(I12,RHO1,RHO2,
     &                  IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                  NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  NAEL,NBEL,
     &                  IJAGRP,IJBGRP,
     &                  SB,CB,C2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  X,NSMOB,NSMST,NSMSX,NSMDX,
     &                  NIA,NIB,NJA,NJB,MXPOBS,IPRNT,NACOB,RHO1S,
     &                  ieaw,n1,n2)
*
* Contributions to density matrix from sigma block (iasm iatp, ibsm ibtp ) and
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
* NGAS : Number of As'es
* IAOC : Occpation of each AS for alpha strings in L
* IBOC : Occpation of each AS for beta  strings in L
* JAOC : Occpation of each AS for alpha strings in R
* JBOC : Occpation of each AS for beta  strings in R
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* MXPNGAS : Largest number of As'es allowed by program
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ieaw=0 Singlet
* ieaw=1 Triplet
* ======
* Output
* ======
* Rho1, RHo2 : Updated density blocks
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
      INTEGER ADSXA(*)
      INTEGER SXSTST(*),STSTSX(*),DXSTST(*),STSTDX(*),SXDXSX(*)
*. Input
      DIMENSION CB(*),SB(*)
*. Output
      DIMENSION RHO1(n1,*),RHO2(n2,*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION  I1(*),XI1S(*),I2(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION C2(*),RHO1S(*),X(*)
      DIMENSION IAOC(*),IBOC(*),JAOC(*),JBOC(*),NOBPTS(*),IOBPTS(*)
      DIMENSION ITSOB(1)
*
      IACTIVE = 0
*
      iUseab=0
      ii=1
      If (ieaw.eq.1) ii=2
      IF(NBEL.GE.1.AND.IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
*  beta contribution to RHO1
* =============================
*
        CALL GSBBD1(RHO1(1,ii),
     &       NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &       NGAS,IBOC,JBOC,
     &       SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,X(1),
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S)
*
* ================================
* beta-beta contribution to RHO2
* ================================
*
        ii=1
        If (ieaw.eq.1) ii=2
        IF(I12.EQ.2.AND.NBEL.GE.2) THEN
          CALL GSBBD2A(RHO2(1,ii),
     &         NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &         NGAS,IBOC,JBOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS)
*
        END IF
      END IF
*
      IF(NAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* =============================
*  alpha contribution to RHO1
* =============================
*
        ii=1
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
        CALL GSBBD1(RHO1(1,ii),
     &                   NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &                   NGAS,IAOC,JAOC,
     &                   SB,CB,
     &                   ADSXA,SXSTST,STSTSX,MXPNGAS,
     &                   NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,X(1),
     &                   NSMOB,NSMST,NSMSX,MXPOBS,RHO1S)
*
* ===================================
*  alpha-alpha contribution to RHO2
* ===================================
*
        ii=1
        IF(I12.EQ.2.AND.NAEL.GE.2) THEN
          CALL GSBBD2A(RHO2(1,ii),
     &         NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &         NGAS,IAOC,JAOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS)
        END IF
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMAT(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
* ===================================
*  alpha-beta contribution to RHO2
* ===================================
*
      ii=1
      If (ieaw.eq.1) ii=3
      IF(I12.EQ.2.AND.NAEL.GE.1.AND.NBEL.GE.1) THEN
*. Routine uses transposed blocks
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
        CALL GSBBD2B(RHO2(1,ii),
     &                    IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                    JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                    IJAGRP,IJBGRP,NGAS,
     &                    IAOC(1),IBOC(1),JAOC(1),JBOC(1),
     &                    SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                    NOBPTS,IOBPTS,MAXK,
     &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                    NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                    SSCR,CSCR,NACOB,NTEST,ieaw)
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMAT(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        Call Unused_integer_array(DXSTST)
        Call Unused_integer_array(STSTDX)
        Call Unused_integer(IPRNT)
      END IF
      END
