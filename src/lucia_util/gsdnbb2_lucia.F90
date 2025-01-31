!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

subroutine GSDNBB2_LUCIA(I12,RHO1,RHO2,RHO2S,RHO2A,IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,NGAS,IAOC,IBOC,JAOC,JBOC,NAEL,NBEL, &
                         IJAGRP,IJBGRP,SB,CB,C2,ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR, &
                         CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSMOB,NSMST,NSMSX,NSMDX,NIA,NIB,NJA,NJB,MXPOBS,IPRNT,NACOB,RHO1S, &
                         SCLFAC,S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IPACK)
! SUBROUTINE GSDNBB2_LUCIA --> 66
!
! Contributions to density matrix from sigma block (iasm iatp, ibsm ibtp) and
! C block (jasm jatp, jbsm, jbtp)
!
! =====
! Input
! =====
!
! IASM,IATP : Symmetry and type of alpha strings in sigma
! IBSM,IBTP : Symmetry and type of beta  strings in sigma
! JASM,JATP : Symmetry and type of alpha strings in C
! JBSM,JBTP : Symmetry and type of beta  strings in C
! NGAS : Number of As'es
! IAOC : Occpation of each AS for alpha strings in L
! IBOC : Occpation of each AS for beta  strings in L
! JAOC : Occpation of each AS for alpha strings in R
! JBOC : Occpation of each AS for beta  strings in R
! NAEL : Number of alpha electrons
! NBEL : Number of  beta electrons
! IJAGRP    : IA and JA belongs to this group of strings
! IJBGRP    : IB and JB belongs to this group of strings
! CB : Input c block
! ADASX : sym of a+, a => sym of a+a
! ADSXA : sym of a+, a+a => sym of a
! SXSTST : Sym of sx,!st> => sym of sx !st>
! STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
!          is nonvanishing by symmetry
! DXSTST : Sym of dx,!st> => sym of dx !st>
! STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
!          is nonvanishing by symmetry
! MXPNGAS : Largest number of As'es allowed by program
! NOBPTS  : Number of orbitals per type and symmetry
! IOBPTS : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! MAXI   : Largest Number of ' spectator strings 'treated simultaneously
! MAXK   : Largest number of inner resolution strings treated at simult.
! IPACK  : Should 2-body densities be packed?
!
! ======
! Output
! ======
! Rho1, RHo2, RHo2s, RHo2a : Updated density blocks
!
! =======
! Scratch
! =======
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
!              largest number of orbital pairs of given symmetries and
!              types.
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! C2 : Must hold largest STT block of sigma or C
!
! XINT : Scratch space for integrals.
!
! Jeppe Olsen, Winter of 1991

implicit real*8(A-H,O-Z)
#include "timers.fh"
integer ADSXA(*), SXSTST(*), STSTSX(*), DXSTST(*), STSTDX(*), SXDXSX(*)
! Input
dimension CB(*), SB(*), NOBPTS(*), IOBPTS(*), X(*)
logical IPACK
! Output
dimension RHO1(*), RHO2(*), RHO2S(*), RHO2A(*)
! Scratch
dimension SSCR(*), CSCR(*)
dimension I1(*), XI1S(*), I2(*), XI2S(*), I3(*), XI3S(*), I4(*), XI4S(*)
dimension C2(*), RHO1S(*)
dimension IAOC(*), JAOC(*), IBOC(*), JBOC(*)
dimension ITSOB(1), IPHGAS(*), SRHO1(*)

NTEST = 0
NTEST = max(NTEST,IPRNT)
if (NTEST >= 200) then
  write(6,*) ' =================='
  write(6,*) ' GSDNBB2 :  R block'
  write(6,*) ' =================='
  call WRTMAT(CB,NJA,NJB,NJA,NJB)
  write(6,*) ' =================='
  write(6,*) ' GSDNBB2 :  L block'
  write(6,*) ' =================='
  call WRTMAT(SB,NIA,NIB,NIA,NIB)

  write(6,*)
  write(6,*) ' Occupation of alpha strings in L'
  call IWRTMA(IAOC,1,NGAS,1,NGAS)
  write(6,*)
  write(6,*) ' Occupation of beta  strings in L'
  call IWRTMA(IBOC,1,NGAS,1,NGAS)
  write(6,*)
  write(6,*) ' Occupation of alpha strings in R'
  call IWRTMA(JAOC,1,NGAS,1,NGAS)
  write(6,*)
  write(6,*) ' Occupation of beta  strings in R'
  call IWRTMA(JBOC,1,NGAS,1,NGAS)

  write(6,*) ' MAXI,MAXK,NSMOB',MAXI,MAXK,NSMOB

  write(6,*) 'SCLFAC =',SCLFAC
end if

if ((IATP == JATP) .and. (IASM == JASM)) then

  ! =========================
  ! beta contribution to RHO1
  ! =========================

  !write(6,*) ' GSBBD1 will be called (beta)'
  IAB = 2
  call TIMING(CPU0,CPU,WALL0,WALL)
  call GSBBD1_LUCIA(RHO1,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NGAS,IBOC,JBOC,SB,CB,ADSXA,SXSTST,STSTSX,MXPNGAS,NOBPTS,IOBPTS, &
                    ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,X,NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,IUSE_PH,IPHGAS,IDOSRHO1, &
                    SRHO1,IAB)
  call TIMING(CPU1,CPU,WALL1,WALL)
  TDENSI(1) = TDENSI(1)+(WALL1-WALL0)

  ! CALL GSBBD1_LUCIA --> 40

  ! ==============================
  ! beta-beta contribution to RHO2
  ! ==============================

  if ((I12 == 2) .and. (NBEL >= 2)) then
    !write(6,*) ' GSBBD2A will be called (beta)'
    call TIMING(CPU0,CPU,WALL0,WALL)
    call GSBBD2A_LUCIA(RHO2,RHO2S,RHO2A,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NGAS,IBOC,JBOC,SB,CB,ADSXA,SXSTST,STSTSX,SXDXSX, &
                       MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,X,NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,IPACK)
    !    GSBBD2A_LUCIA(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS, &
    !                  NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,X,NSMOB,NSMST,NSMSX,MXPOBS)
    call TIMING(CPU1,CPU,WALL1,WALL)
    TDENSI(2) = TDENSI(2)+(WALL1-WALL0)

    ! CALL GSBBD2A_LUCIA --> 37

    !write(6,*) ' GSBBD2A was called'

  end if
end if

if ((IBTP == JBTP) .and. (IBSM == JBSM)) then

  ! ==========================
  ! alpha contribution to RHO1
  ! ==========================

  call TRPMT3(CB,NJA,NJB,C2)
  call COPVEC(C2,CB,NJA*NJB)
  call TRPMT3(SB,NIA,NIB,C2)
  call COPVEC(C2,SB,NIA*NIB)
  !write(6,*) ' GSBBD1 will be called (alpha)'
  IAB = 1
  call TIMING(CPU0,CPU,WALL0,WALL)
  call GSBBD1_LUCIA(RHO1,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,NGAS,IAOC,JAOC,SB,CB,ADSXA,SXSTST,STSTSX,MXPNGAS,NOBPTS,IOBPTS, &
                    ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,X,NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,IUSE_PH,IPHGAS,IDOSRHO1, &
                    SRHO1,IAB)
  call TIMING(CPU1,CPU,WALL1,WALL)
  TDENSI(1) = TDENSI(1)+(WALL1-WALL0)

  ! CALL GSBBD1_LUCIA --> 40

  !write(6,*) ' GSBBD1 was called'
  if ((I12 == 2) .and. (NAEL >= 2)) then

    ! ================================
    ! alpha-alpha contribution to RHO2
    ! ================================

    !write(6,*) ' GSBBD2A will be called (alpha)'
    call TIMING(CPU0,CPU,WALL0,WALL)
    call GSBBD2A_LUCIA(RHO2,RHO2S,RHO2A,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,NGAS,IAOC,JAOC,SB,CB,ADSXA,SXSTST,STSTSX,SXDXSX, &
                       MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,X,NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,IPACK)
    call TIMING(CPU1,CPU,WALL1,WALL)
    TDENSI(2) = TDENSI(2)+(WALL1-WALL0)

    ! CALL GSBBD2A_LUCIA --> 37

    !write(6,*) ' GSBBD2A was called'
  end if
  call TRPMT3(CB,NJB,NJA,C2)
  call COPVEC(C2,CB,NJA*NJB)
  call TRPMAT(SB,NIB,NIA,C2)
  call COPVEC(C2,SB,NIB*NIA)
end if
!
! ===============================
! alpha-beta contribution to RHO2
! ===============================
!
if ((I12 == 2) .and. (NAEL >= 1) .and. (NBEL >= 1)) then
  ! Routine uses transposed blocks
  call TRPMT3(CB,NJA,NJB,C2)
  call COPVEC(C2,CB,NJA*NJB)
  call TRPMT3(SB,NIA,NIB,C2)
  call COPVEC(C2,SB,NIA*NIB)
  !write(6,*) ' GSBBD2B will be called'
  IUSEAB = 0
  call TIMING(CPU0,CPU,WALL0,WALL)
  call GSBBD2B_LUCIA(RHO2,RHO2S,RHO2A,IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC, &
                     JBOC,SB,CB,ADSXA,STSTSX,MXPNGAS,NOBPTS,IOBPTS,MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSMOB,NSMST,NSMSX,NSMDX, &
                     MXPOBS,IUSEAB,SSCR,CSCR,NACOB,NTEST,SCLFAC,S2_TERM1,IPACK)
  !    GSBBD2B_LUCIA(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,SB,CB, &
  !                  ADSXA,STSTSX,MXPNGAS,NOBPTS,IOBPTS,MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS, &
  !                  IUSEAB,CJRES,SIRES,NORB,NTEST)
  call TIMING(CPU1,CPU,WALL1,WALL)
  TDENSI(3) = TDENSI(3)+(WALL1-WALL0)

  ! CALL GSBBD2B_LUCIA --> 52

  !write(6,*) ' GSBBD2B was called'

  call TRPMT3(CB,NJB,NJA,C2)
  call COPVEC(C2,CB,NJA*NJB)
  call TRPMAT(SB,NIB,NIA,C2)
  call COPVEC(C2,SB,NIB*NIA)
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(DXSTST)
  call Unused_integer_array(STSTDX)
end if

end subroutine GSDNBB2_LUCIA
