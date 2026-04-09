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

!#define _DEBUGPRINT_
subroutine GSDNBB2(I12,RHO1,RHO2,RHO2S,RHO2A,IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,NGAS,IAOC,IBOC,JAOC,JBOC,NAEL,NBEL,IJAGRP, &
                   IJBGRP,SB,CB,C2,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSMOB,NIA,NIB,NJA, &
                   NJB,NACOB,RHO1S,SCLFAC,S2_TERM1,IPHGAS,IDOSRHO1,SRHO1,IPACK)
! SUBROUTINE GSDNBB2 --> 66
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
! MXPNGAS : Largest number of As'es allowed by program
! NOBPTS  : Number of orbitals per type and symmetry
! IOBPTS : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! MAXI   : Largest Number of "spectator strings" treated simultaneously
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

use lucia_data, only: TDENSI
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: I12, IASM, IATP, IBSM, IBTP, JASM, JATP, JBSM, JBTP, NGAS, IAOC(*), IBOC(*), JAOC(*), JBOC(*), &
                                 NAEL, NBEL, IJAGRP, IJBGRP, MXPNGAS, NOBPTS(*), IOBPTS(*), MAXI, MAXK, NSMOB, NIA, NIB, NJA, NJB, &
                                 NACOB, IPHGAS(*), IDOSRHO1
real(kind=wp), intent(inout) :: RHO1(*), RHO2(*), RHO2S(*), RHO2A(*), SB(*), CB(*), XI1S(*), XI2S(*), RHO1S(*), S2_TERM1, SRHO1(*)
real(kind=wp), intent(_OUT_) :: C2(*), SSCR(*), CSCR(*), XI3S(*), XI4S(*), X(*)
integer(kind=iwp), intent(inout) :: I1(*), I2(*)
integer(kind=iwp), intent(_OUT_) :: I3(*), I4(*)
real(kind=wp), intent(_OUT_) :: SCLFAC
logical(kind=iwp), intent(in) :: IPACK
real(kind=wp) :: CPU, CPU0, CPU1, WALL, WALL0, WALL1
integer(kind=iwp) :: IAB, IUSEAB

#ifdef _DEBUGPRINT_
write(u6,*) ' =================='
write(u6,*) ' GSDNBB2 :  R block'
write(u6,*) ' =================='
call WRTMAT(CB,NJA,NJB,NJA,NJB)
write(u6,*) ' =================='
write(u6,*) ' GSDNBB2 :  L block'
write(u6,*) ' =================='
call WRTMAT(SB,NIA,NIB,NIA,NIB)

write(u6,*)
write(u6,*) ' Occupation of alpha strings in L'
call IWRTMA(IAOC,1,NGAS,1,NGAS)
write(u6,*)
write(u6,*) ' Occupation of beta  strings in L'
call IWRTMA(IBOC,1,NGAS,1,NGAS)
write(u6,*)
write(u6,*) ' Occupation of alpha strings in R'
call IWRTMA(JAOC,1,NGAS,1,NGAS)
write(u6,*)
write(u6,*) ' Occupation of beta  strings in R'
call IWRTMA(JBOC,1,NGAS,1,NGAS)

write(u6,*) ' MAXI,MAXK,NSMOB',MAXI,MAXK,NSMOB

write(u6,*) 'SCLFAC =',SCLFAC
#endif

if ((IATP == JATP) .and. (IASM == JASM)) then

  ! =========================
  ! beta contribution to RHO1
  ! =========================

  !write(u6,*) ' GSBBD1 will be called (beta)'
  IAB = 2
  call TIMING(CPU0,CPU,WALL0,WALL)
  call GSBBD1(RHO1,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NGAS,IBOC,JBOC,SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2, &
              XI2S,NSMOB,RHO1S,SCLFAC,IPHGAS,IDOSRHO1,SRHO1,IAB)
  call TIMING(CPU1,CPU,WALL1,WALL)
  TDENSI(1) = TDENSI(1)+(WALL1-WALL0)

  ! CALL GSBBD1 --> 40

  ! ==============================
  ! beta-beta contribution to RHO2
  ! ==============================

  if ((I12 == 2) .and. (NBEL >= 2)) then
    !write(u6,*) ' GSBBD2A will be called (beta)'
    call TIMING(CPU0,CPU,WALL0,WALL)
    call GSBBD2A(RHO2,RHO2S,RHO2A,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NGAS,IBOC,JBOC,SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR, &
                 CSCR,I1,XI1S,X,NSMOB,SCLFAC,IPACK)
    !    GSBBD2A(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1, &
    !            XI1S,X,NSMOB)
    call TIMING(CPU1,CPU,WALL1,WALL)
    TDENSI(2) = TDENSI(2)+(WALL1-WALL0)

    ! CALL GSBBD2A --> 37

    !write(u6,*) ' GSBBD2A was called'

  end if
end if

if ((IBTP == JBTP) .and. (IBSM == JBSM)) then

  ! ==========================
  ! alpha contribution to RHO1
  ! ==========================

  call TRPMT3(CB,NJA,NJB,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRPMT3(SB,NIA,NIB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
  !write(u6,*) ' GSBBD1 will be called (alpha)'
  IAB = 1
  call TIMING(CPU0,CPU,WALL0,WALL)
  call GSBBD1(RHO1,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,NGAS,IAOC,JAOC,SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2, &
              XI2S,NSMOB,RHO1S,SCLFAC,IPHGAS,IDOSRHO1,SRHO1,IAB)
  call TIMING(CPU1,CPU,WALL1,WALL)
  TDENSI(1) = TDENSI(1)+(WALL1-WALL0)

  ! CALL GSBBD1 --> 40

  !write(u6,*) ' GSBBD1 was called'
  if ((I12 == 2) .and. (NAEL >= 2)) then

    ! ================================
    ! alpha-alpha contribution to RHO2
    ! ================================

    !write(u6,*) ' GSBBD2A will be called (alpha)'
    call TIMING(CPU0,CPU,WALL0,WALL)
    call GSBBD2A(RHO2,RHO2S,RHO2A,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,NGAS,IAOC,JAOC,SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR, &
                 CSCR,I1,XI1S,X,NSMOB,SCLFAC,IPACK)
    call TIMING(CPU1,CPU,WALL1,WALL)
    TDENSI(2) = TDENSI(2)+(WALL1-WALL0)

    ! CALL GSBBD2A --> 37

    !write(u6,*) ' GSBBD2A was called'
  end if
  call TRPMT3(CB,NJB,NJA,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRNSPS(NIB,NIA,SB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
end if

! ===============================
! alpha-beta contribution to RHO2
! ===============================

if ((I12 == 2) .and. (NAEL >= 1) .and. (NBEL >= 1)) then
  ! Routine uses transposed blocks
  call TRPMT3(CB,NJA,NJB,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRPMT3(SB,NIA,NIB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
  !write(u6,*) ' GSBBD2B will be called'
  IUSEAB = 0
  call TIMING(CPU0,CPU,WALL0,WALL)
  call GSBBD2B(RHO2,RHO2S,RHO2A,IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,SB, &
               CB,MXPNGAS,NOBPTS,IOBPTS,MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSMOB,IUSEAB,SSCR,CSCR,NACOB,SCLFAC,S2_TERM1,IPACK)
  !    GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,SB,CB,MXPNGAS, &
  !            NOBPTS,IOBPTS,MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSMOB,IUSEAB,CJRES,SIRES,NORB)
  call TIMING(CPU1,CPU,WALL1,WALL)
  TDENSI(3) = TDENSI(3)+(WALL1-WALL0)

  ! CALL GSBBD2B --> 52

  !write(u6,*) ' GSBBD2B was called'

  call TRPMT3(CB,NJB,NJA,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRNSPS(NIB,NIA,SB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
end if

end subroutine GSDNBB2
