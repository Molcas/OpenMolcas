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

subroutine GSDNBB2_MCLR(I12,RHO1,RHO2,IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,NGAS,IAOC,IBOC,JAOC,JBOC,NAEL,NBEL,IJAGRP,IJBGRP,SB, &
                        CB,C2,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSM,NIA,NIB,NJA,NJB, &
                        NACOB,RHO1S,ieaw,n1,n2)
! Contributions to density matrix from sigma block (iasm iatp, ibsm ibtp) and
! C block (jasm jatp, jbsm, jbtp)
!
! =====
! Input
! =====
! IASM,IATP : Symmetry and type of alpha strings in sigma
! IBSM,IBTP : Symmetry and type of beta  strings in sigma
! JASM,JATP : Symmetry and type of alpha strings in C
! JBSM,JBTP : Symmetry and type of beta  strings in C
! NGAS      : Number of As'es
! IAOC      : Occpation of each AS for alpha strings in L
! IBOC      : Occpation of each AS for beta  strings in L
! JAOC      : Occpation of each AS for alpha strings in R
! JBOC      : Occpation of each AS for beta  strings in R
! NAEL      : Number of alpha electrons
! NBEL      : Number of  beta electrons
! IJAGRP    : IA and JA belongs to this group of strings
! IJBGRP    : IB and JB belongs to this group of strings
! CB        : Input c block
! MXPNGAS   : Largest number of As'es allowed by program
! NOBPTS    : Number of orbitals per type and symmetry
! IOBPTS    : base for orbitals of given type and symmetry
! IBORB     : Orbitals of given type and symmetry
! MAXI      : Largest Number of "spectator strings" treated simultaneously
! MAXK      : Largest number of inner resolution strings treated at simult.
!
! ieaw=0 Singlet
! ieaw=1 Triplet

! ======
! Output
! ======
! Rho1, RHo2 : Updated density blocks

! =======
! Scratch
! =======
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the largest
!              number of orbital pairs of given symmetries and types.
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! C2         : Must hold largest STT block of sigma or C
! XINT       : Scratch space for integrals.
!
! Jeppe Olsen, Winter of 1991

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: I12, IASM, IATP, IBSM, IBTP, JASM, JATP, JBSM, JBTP, IAOC(*), IBOC(*), JAOC(*), JBOC(*), NAEL, &
                                 NBEL, IJAGRP, IJBGRP, MXPNGAS, NOBPTS(*), IOBPTS(*), MAXI, MAXK, NSM, NIA, NIB, NJA, NJB, NACOB, &
                                 ieaw, n1, n2
real(kind=wp), intent(inout) :: RHO1(n1,*), RHO2(n2,*), SB(*), CB(*)
integer(kind=iwp), intent(inout) :: NGAS
real(kind=wp), intent(_OUT_) :: C2(*), SSCR(*), CSCR(*), XI1S(*), XI2S(*), XI3S(*), XI4S(*), X(*), RHO1S(*)
integer(kind=iwp), intent(_OUT_) :: I1(*), I2(*), I3(*), I4(*)
integer(kind=iwp) :: ii, iUseab

iUseab = 0
ii = 1
if (ieaw == 1) ii = 2
if ((NBEL >= 1) .and. (IATP == JATP) .and. (JASM == IASM)) then

  ! ===========================
  !  beta contribution to RHO1
  ! ===========================

  call GSBBD1_MCLR(RHO1(:,ii),NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NGAS,IBOC,JBOC,SB,CB,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S, &
                   I2,XI2S,NSM,RHO1S)

  ! ================================
  !  beta-beta contribution to RHO2
  ! ================================

  ii = 1
  if (ieaw == 1) ii = 2
  if ((I12 == 2) .and. (NBEL >= 2)) &
    call GSBBD2A_MCLR(RHO2(:,ii),NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NGAS,IBOC,JBOC,SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR, &
                      CSCR,I1,XI1S,X,NSM)

end if

if ((NAEL >= 1) .and. (IBTP == JBTP) .and. (IBSM == JBSM)) then

  ! ============================
  !  alpha contribution to RHO1
  ! ============================

  ii = 1
  call TRPMT3(CB,NJA,NJB,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRPMT3(SB,NIA,NIB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
  call GSBBD1_MCLR(RHO1(:,ii),NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,NGAS,IAOC,JAOC,SB,CB,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S, &
                   I2,XI2S,NSM,RHO1S)

  ! ==================================
  !  alpha-alpha contribution to RHO2
  ! ==================================

  ii = 1
  if ((I12 == 2) .and. (NAEL >= 2)) &
    call GSBBD2A_MCLR(RHO2(:,ii),NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,NGAS,IAOC,JAOC,SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR, &
                      CSCR,I1,XI1S,X,NSM)
  call TRPMT3(CB,NJB,NJA,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRNSPS(NIB,NIA,SB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
end if

! =================================
!  alpha-beta contribution to RHO2
! =================================

ii = 1
if (ieaw == 1) ii = 3
if ((I12 == 2) .and. (NAEL >= 1) .and. (NBEL >= 1)) then
  ! Routine uses transposed blocks
  call TRPMT3(CB,NJA,NJB,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRPMT3(SB,NIA,NIB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
  call GSBBD2B_MCLR(RHO2(:,ii),IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,SB, &
                    CB,NOBPTS,IOBPTS,MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSM,IUSEAB,SSCR,CSCR,NACOB,ieaw)
  call TRPMT3(CB,NJB,NJA,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRNSPS(NIB,NIA,SB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
end if

return

end subroutine GSDNBB2_MCLR
