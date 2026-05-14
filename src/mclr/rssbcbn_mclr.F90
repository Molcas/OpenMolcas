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

subroutine RSSBCBN_MCLR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,IAEL1,IAEL3,IBEL1,IBEL3,JAEL1,JAEL3,JBEL1,JBEL3,NAEL,NBEL,IJAGRP, &
                        IJBGRP,SB,CB,IDOH2,NTSOB,IBTSOB,ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,C2,NSM,NIA, &
                        NIB,NJA,NJB,IST,CJRES,SIRES,NOPART,TimeDep)
! Contributions to sigma block (iasm iatp, ibsm ibtp) from
! C block (jasm jatp, jbsm, jbtp)
!
! =====
! Input
! =====
! IASM,IATP   : Symmetry and type of alpha strings in sigma
! IBSM,IBTP   : Symmetry and type of beta  strings in sigma
! JASM,JATP   : Symmetry and type of alpha strings in C
! JBSM,JBTP   : Symmetry and type of beta  strings in C
! IAEL1,IAEL3 : Number of elecs in RAS1(RAS3) for alpha strings in sigma
! IBEL1,IBEL3 : Number of elecs in RAS1(RAS3) for  beta strings in sigma
! JAEL1,JAEL3 : Number of elecs in RAS1(RAS3) for alpha strings in C
! JBEL1,JBEL3 : Number of elecs in RAS1(RAS3) for  beta strings in C
! NAEL        : Number of alpha electrons
! NBEL        : Number of  beta electrons
! IJAGRP      : IA and JA belongs to this group of strings
! IJBGRP      : IB and JB belongs to this group of strings
! CB          : Input c block
! IDOH2       : = 0 => no two electron operator
! IDOH2       : = 1 =>    two electron operator
! NTSOB       : Number of orbitals per type and symmetry
! IBTSOB      : base for orbitals of given type and symmetry
! IBORB       : Orbitals of given type and symmetry
! MAXI        : Largest Number of "spectator strings" treated simultaneously
! MAXK        : Largest number of inner resolution strings treated at simult.
! IST, IDOH2  : See RASSG3 input description
!
! ======
! Output
! ======
! SB : fresh sigma block
!
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

use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IASM, IATP, IBSM, IBTP, JASM, JATP, JBSM, JBTP, IAEL1, IAEL3, IBEL1, IBEL3, JAEL1, JAEL3, JBEL1, &
                                 JBEL3, NAEL, NBEL, IJAGRP, IJBGRP, IDOH2, NTSOB(*), IBTSOB(*), ITSOB(*), MAXI, MAXK, NSM, NIA, &
                                 NIB, NJA, NJB, IST, NOPART
real(kind=wp), intent(inout) :: SB(*), CB(*)
real(kind=wp), intent(_OUT_) :: SSCR(*), CSCR(*), XI1S(MAXK,*), XI2S(MAXK,*), XI3S(MAXK,*), XI4S(MAXK,*), XINT(*), C2(*), &
                                CJRES(*), SIRES(*)
integer(kind=iwp), intent(_OUT_) :: I1(MAXK,*), I2(MAXK,*), I3(MAXK,*), I4(MAXK,*)
logical(kind=iwp), intent(in) :: TimeDep
integer(kind=iwp) :: ieaw, IFACTOR, IIITRNS, JJJTRNS
real(kind=wp) :: SGN

! ============================
! Sigma beta beta contribution
! ============================

! Sigma aa(IA,IB) = sum(i > k,j > l)<IB!Eb(ij)Eb(kl)!JB>
!                 * ((ij!kl)-(il!kj)) C(IA,JB)
!                 + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)

!write(u6,*) 'I am in rssbcbn'
if ((IATP == JATP) .and. (JASM == IASM)) then

  ! One electron part

  if (IST == 1) then
    SGN = One
  else
    SGN = -One
  end if
  if (NBEL >= 1) &
    call RSBB1E_MCLR(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,IBEL1,IBEL3,JBEL1,JBEL3,SB,CB,NTSOB,IBTSOB,ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S, &
                     XINT,NSM,SGN)

  ! Two electron part

  if ((IDOH2 /= 0) .and. (NBEL >= 2)) then
    !write(u6,*) 'Timedep in rssbcbn',TimeDep
    ieaw = 0
    call RSBB2A_MCLR(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,IBEL1,IBEL3,JBEL1,JBEL3,SB,CB,NTSOB,IBTSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,XINT, &
                     NSM,SGN,NOPART,TimeDep,ieaw)
  end if
end if

!====================================
! Mixed alpha-beta double excitations
!====================================

if ((IDOH2 /= 0) .and. (NAEL >= 1) .and. (NBEL >= 1)) then

  ieaw = 0
  if (ist == 2) ieaw = 1
  call TRNSPS(NJA,NJB,CB,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRNSPS(NIA,NIB,SB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
  IIITRNS = 1
  if ((IIITRNS == 1) .and. (NIB > NIA) .and. (NJB > NJA)) then
    JJJTRNS = 1
  else
    JJJTRNS = 0
  end if
  if ((JJJTRNS == 1) .and. (IST == 2)) then
    IFACTOR = -1
  else
    IFACTOR = 1
  end if
  if (JJJTRNS == 0) then
    call RSBB2BN_MCLR(IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IJAGRP,IJBGRP,IAEL1,IAEL3,JAEL1,JAEL3,IBEL1,IBEL3, &
                      JBEL1,JBEL3,SB,CB,NTSOB,IBTSOB,MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,NSM,CJRES,SIRES,IFACTOR,ieaw,TimeDep)
  else if (JJJTRNS == 1) then
    call TRNSPS(NIB,NIA,SB,C2)
    SB(1:NIA*NIB) = C2(1:NIA*NIB)
    call TRNSPS(NJB,NJA,CB,C2)
    CB(1:NJA*NJB) = C2(1:NJA*NJB)

    call RSBB2BN_MCLR(IBSM,IBTP,IASM,IATP,NIB,NIA,JBSM,JBTP,JASM,JATP,NJB,NJA,IJBGRP,IJAGRP,IBEL1,IBEL3,JBEL1,JBEL3,IAEL1,IAEL3, &
                      JAEL1,JAEL3,SB,CB,NTSOB,IBTSOB,MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,NSM,CJRES,SIRES,IFACTOR,ieaw,TimeDep)
    call TRNSPS(NIA,NIB,SB,C2)
    SB(1:NIA*NIB) = C2(1:NIA*NIB)
    call TRNSPS(NJA,NJB,CB,C2)
    CB(1:NJA*NJB) = C2(1:NJA*NJB)
  end if
  ! Restore order!
  call TRNSPS(NJB,NJA,CB,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRNSPS(NIB,NIA,SB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
end if

! ========================
! Sigma alpha contribution
! ========================

! Transpose for alpha excitations

if ((NAEL >= 1) .and. (IBTP == JBTP) .and. (IBSM == JBSM)) then
  call TRNSPS(NJA,NJB,CB,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
  call TRNSPS(NIA,NIB,SB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)

  ! alpha single excitation

  SGN = One
  call RSBB1E_MCLR(IASM,IATP,JASM,JATP,IJAGRP,NIB,IAEL1,IAEL3,JAEL1,JAEL3,SB,CB,NTSOB,IBTSOB,ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S, &
                   XINT,NSM,SGN)

  ! alpha double excitation

  if ((NAEL >= 2) .and. (IDOH2 /= 0)) &
    call RSBB2A_MCLR(IASM,IATP,JASM,JATP,IJAGRP,NIB,IAEL1,IAEL3,JAEL1,JAEL3,SB,CB,NTSOB,IBTSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,XINT, &
                     NSM,SGN,NOPART,TimeDep,ieaw)

  ! Restore order!
  call TRNSPS(NIB,NIA,SB,C2)
  SB(1:NIA*NIB) = C2(1:NIA*NIB)
  call TRNSPS(NJB,NJA,CB,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
end if

return

end subroutine RSSBCBN_MCLR
