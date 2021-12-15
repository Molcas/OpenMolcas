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
      SUBROUTINE GSDNBB2_LUCIA(    I12,   RHO1,   RHO2,  RHO2S,  RHO2A,
     &                            IASM,   IATP,   IBSM,   IBTP,   JASM,
     &                            JATP,   JBSM,   JBTP,   NGAS,   IAOC,
     &                            IBOC,   JAOC,   JBOC,   NAEL,   NBEL,
     &                          IJAGRP, IJBGRP,     SB,     CB,     C2,
*
     &                           ADSXA, SXSTST, STSTSX, DXSTST, STSTDX,
     &                          SXDXSX,MXPNGAS, NOBPTS, IOBPTS,   MAXI,
     &                            MAXK,   SSCR,   CSCR,     I1,   XI1S,
     &                              I2,   XI2S,     I3,   XI3S,     I4,
     &                            XI4S,      X,  NSMOB,  NSMST,  NSMSX,
*
     &                           NSMDX,    NIA,    NIB,    NJA,    NJB,
     &                          MXPOBS,  IPRNT,  NACOB,  RHO1S, SCLFAC,
     &                         S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,
     &                           IPACK)
*
* SUBROUTINE GSDNBB2_LUCIA --> 66
*
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
* IPACK  : Should 2-body densities be packed?
*
* ======
* Output
* ======
* Rho1, RHo2, RHo2s, RHo2a : Updated density blocks
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
      INTEGER ADSXA(*),SXSTST(*),STSTSX(*),DXSTST(*),STSTDX(*),SXDXSX(*)
*. Input
      DIMENSION CB(*),SB(*),NOBPTS(*),IOBPTS(*),X(*)
      LOGICAL   IPACK
*. Output
      DIMENSION RHO1(*),RHO2(*),RHO2S(*),RHO2A(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION  I1(*),XI1S(*),I2(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION C2(*),RHO1S(*)

      DIMENSION IAOC(*),JAOC(*),IBOC(*),JBOC(*)
      DIMENSION ITSOB(1),IPHGAS(*),SRHO1(*)
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      NTESTO= NTEST
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' =================='
        WRITE(6,*) ' GSDNBB2 :  R block '
        WRITE(6,*) ' ==================='
        CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' GSDNBB2 :  L block '
        WRITE(6,*) ' ==================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
*
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in L '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in L '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in R '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in R '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' MAXI,MAXK,NSMOB',MAXI,MAXK,NSMOB
*
        WRITE(6,*) 'SCLFAC =',SCLFAC
      END IF
      IACTIVE = 0
*
      IF(IATP.EQ.JATP.AND.IASM.EQ.JASM) THEN
*
* =============================
*  beta contribution to RHO1
* =============================
*
C?      WRITE(6,*) ' GSBBD1 will be called (beta)'
        IAB = 2
        CALL TIMING(CPU0,CPU,WALL0,WALL)
        CALL GSBBD1_LUCIA(    RHO1,   NACOB,    IBSM,    IBTP,    JBSM,
     &                        JBTP,  IJBGRP,     NIA,    NGAS,    IBOC,
     &                        JBOC,      SB,      CB,   ADSXA,  SXSTST,
     &                      STSTSX, MXPNGAS,  NOBPTS,  IOBPTS,   ITSOB,
     &                        MAXI,    MAXK,    SSCR,    CSCR,      I1,
*
     &                        XI1S,      I2,    XI2S,       X,   NSMOB,
     &                       NSMST,   NSMSX,  MXPOBS,   RHO1S,  SCLFAC,
     &                     IUSE_PH,  IPHGAS,IDOSRHO1,   SRHO1,     IAB)
        CALL TIMING(CPU1,CPU,WALL1,WALL)
        TDENSI(1)=TDENSI(1)+(WALL1-WALL0)
*
* CALL GSBBD1_LUCIA --> 40
*
C?      WRITE(6,*) ' GSBBD1 was called '
C?      WRITE(6,*) ' Memory check '
C?      CALL MEMCHK
*
* ================================
* beta-beta contribution to RHO2
* ================================
*
        IF(I12.EQ.2.AND.NBEL.GE.2) THEN
C?        WRITE(6,*) ' GSBBD2A will be called (beta)'
          CALL TIMING(CPU0,CPU,WALL0,WALL)
          CALL GSBBD2A_LUCIA(   RHO2,  RHO2S,  RHO2A,  NACOB,   IBSM,
     &                          IBTP,   JBSM,   JBTP, IJBGRP,    NIA,
     &                          NGAS,   IBOC,   JBOC,     SB,     CB,
     &                         ADSXA, SXSTST, STSTSX, SXDXSX,MXPNGAS,
     &                        NOBPTS, IOBPTS,   MAXI,   MAXK,   SSCR,
*
     &                          CSCR,     I1,   XI1S,     I2,   XI2S,
     &                             X,  NSMOB,  NSMST,  NSMSX, MXPOBS,
     &                        SCLFAC,  IPACK)
          CALL TIMING(CPU1,CPU,WALL1,WALL)
          TDENSI(2)=TDENSI(2)+(WALL1-WALL0)
*
* CALL GSBBD2A_LUCIA --> 37
*
C?        WRITE(6,*) ' GSBBD2A was called '
*
C              GSBBD2A_LUCIA(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
C    &         NGAS,ISEL,ICEL,SB,CB,
C    &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
C    &         NOBPTS,IOBPTS,MAXI,MAXK,
C    &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
C    &         NSMOB,NSMST,NSMSX,MXPOBS)
        END IF
      END IF
*
      IF(IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* =============================
*  alpha contribution to RHO1
* =============================
*
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
C?      WRITE(6,*) ' GSBBD1 will be called (alpha)'
        IAB = 1
        CALL TIMING(CPU0,CPU,WALL0,WALL)
        CALL GSBBD1_LUCIA(    RHO1,   NACOB,    IASM,    IATP,    JASM,
     &                        JATP,  IJAGRP,     NIB,    NGAS,    IAOC,
     &                        JAOC,      SB,      CB,   ADSXA,  SXSTST,
     &                      STSTSX, MXPNGAS,  NOBPTS,  IOBPTS,   ITSOB,
     &                        MAXI,    MAXK,    SSCR,    CSCR,      I1,
*
     &                        XI1S,      I2,    XI2S,       X,   NSMOB,
     &                       NSMST,   NSMSX,  MXPOBS,   RHO1S,  SCLFAC,
     &                     IUSE_PH,  IPHGAS,IDOSRHO1,   SRHO1,     IAB)
        CALL TIMING(CPU1,CPU,WALL1,WALL)
        TDENSI(1)=TDENSI(1)+(WALL1-WALL0)
*
* CALL GSBBD1_LUCIA --> 40
*
C?      WRITE(6,*) ' GSBBD1 was called '
        IF(I12.EQ.2.AND.NAEL.GE.2) THEN
*
* ===================================
*  alpha-alpha contribution to RHO2
* ===================================
*
C?        WRITE(6,*) ' GSBBD2A will be called (alpha)'
          CALL TIMING(CPU0,CPU,WALL0,WALL)
          CALL GSBBD2A_LUCIA(   RHO2,  RHO2S,  RHO2A,  NACOB,   IASM,
     &                          IATP,   JASM,   JATP, IJAGRP,    NIB,
     &                          NGAS,   IAOC,   JAOC,     SB,     CB,
     &                         ADSXA, SXSTST, STSTSX, SXDXSX,MXPNGAS,
     &                        NOBPTS, IOBPTS,   MAXI,   MAXK,   SSCR,
*
     &                          CSCR,     I1,   XI1S,     I2,   XI2S,
     &                             X,  NSMOB,  NSMST,  NSMSX, MXPOBS,
     &                        SCLFAC,  IPACK)
          CALL TIMING(CPU1,CPU,WALL1,WALL)
          TDENSI(2)=TDENSI(2)+(WALL1-WALL0)
*
* CALL GSBBD2A_LUCIA --> 37
*
C?        WRITE(6,*) ' GSBBD2A was called '
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
      IF(I12.EQ.2.AND.NAEL.GE.1.AND.NBEL.GE.1) THEN
*. Routine uses transposed blocks
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
C?      WRITE(6,*) ' GSBBD2B will be called '
        IUSEAB = 0
        CALL TIMING(CPU0,CPU,WALL0,WALL)
        CALL GSBBD2B_LUCIA(    RHO2,   RHO2S,   RHO2A,    IASM,    IATP,
     &                         IBSM,    IBTP,     NIA,     NIB,    JASM,
     &                         JATP,    JBSM,    JBTP,     NJA,     NJB,
     &                       IJAGRP,  IJBGRP,    NGAS,    IAOC,    IBOC,
     &                         JAOC,    JBOC,      SB,      CB,   ADSXA,
*
     &                       STSTSX, MXPNGAS,  NOBPTS,  IOBPTS,    MAXK,
     &                           I1,    XI1S,      I2,    XI2S,      I3,
     &                         XI3S,      I4,    XI4S,       X,   NSMOB,
     &                        NSMST,   NSMSX,   NSMDX,  MXPOBS,  IUSEAB,
     &                         SSCR,    CSCR,   NACOB,   NTEST,  SCLFAC,
*
     &                     S2_TERM1,   IPACK)
        CALL TIMING(CPU1,CPU,WALL1,WALL)
        TDENSI(3)=TDENSI(3)+(WALL1-WALL0)
*
* CALL GSBBD2B_LUCIA --> 52
*
C?      WRITE(6,*) ' GSBBD2B was called '
C    &
C     GSBBD2B_LUCIA(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
C    &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
C    &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
C    &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
C    &                  NOBPTS,IOBPTS,MAXK,
C    &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
C    &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
C    &                  CJRES,SIRES,NORB,NTEST)
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMAT(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(DXSTST)
        CALL Unused_integer_array(STSTDX)
      END IF
      END
