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
! Copyright (C) 1993, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PRMBLK(IDC,ISGV,IASM,IBSM,IATP,IBTP,PS,PL,JATP,JBTP,JASM,JBSM,ISGN,ITRP,NPERM)
! A block of CI coefficients defined by by IATP,IASM,IBTP,IBSM is given
!
! Obtain the number of other blocks that can be obtained by spin
! and relection symmetry.
!
! Jeppe Olsen, July 1993
!
! =====
! Output
! =====
! JATP(I),JASM(I),JBTP(I),JBSM(I) indices for Block I
! NPERM : Number of blocks  that can be obtained
! ITRP(I) = 1 => block should     be transposed
!         = 0 => block should not be transposed
! ISGN   : Sign to multiply previous block with to getnew sign
!
! There are four types of permutations
!
!    operation   *      JASM  *      JBSM  * JATP * JBTP * Iperm * Sign *
!   *********************************************************************
!   * Identity   *      IASM  *      IBSM  * IATP * IBTP *   0   * 1    *
!   * Ml         * ISGV(IASM) * ISGV(IBSM) * IATP * IBTP *   0   * PL   *
!   * Ms         *      IBSM  *      IASM  * IBTP * IATP *   1   * PS   *
!   * Ms+Ml      * ISGV(IBSM) * ISGV(IASM) * IBTP * IATP *   1   * PS PL*
!   *********************************************************************

use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: IDC, ISGV(*), IASM, IBSM, IATP, IBTP
integer(kind=iwp), intent(out) :: JATP(4), JBTP(4), JASM(4), JBSM(4), ISGN(4), ITRP(4), NPERM
real(kind=wp), intent(in) :: PS, PL
integer(kind=iwp) :: IPERM, ISET, KASM, KATP, KBSM, KBTP, KSIGN, KTRP, LPERM, LSIGN, LTRP

! To eliminate some compiler warnings
KASM = 0
KBSM = 0
KATP = 0
KBTP = 0
KSIGN = 0
KTRP = 0
LSIGN = 0
LTRP = 0

NPERM = 0
do IPERM=1,4
  ISET = 0
  if (IPERM == 1) then

    ! Identity operation

    KASM = IASM
    KBSM = IBSM
    KATP = IATP
    KBTP = IBTP
    KSIGN = 1
    KTRP = 0
    ISET = 1
  else if ((IPERM == 2) .and. ((IDC == 3) .or. (IDC == 4))) then

    ! Ml reflection

    KASM = ISGV(IASM)
    KBSM = ISGV(IBSM)
    KATP = IATP
    KBTP = IBTP
    if (PL == One) then
      KSIGN = 1
    else if (PL == -One) then
      KSIGN = -1
    end if
    KTRP = 0
    ISET = 1
  else if ((IPERM == 3) .and. ((IDC == 2) .or. (IDC == 4))) then

    ! Ms reflection

    KASM = IBSM
    KBSM = IASM
    KATP = IBTP
    KBTP = IATP
    if (PS == One) then
      KSIGN = 1
    else if (PS == -One) then
      KSIGN = -1
    end if
    KTRP = 1
    ISET = 1
  else if ((IPERM == 4) .and. (IDC == 4)) then

    ! Ms Ml  reflection

    KASM = ISGV(IBSM)
    KBSM = ISGV(IASM)
    KATP = IBTP
    KBTP = IATP
    if (PS*PL == One) then
      KSIGN = 1
    else if (PS == -One) then
      KSIGN = -1
    end if
    KTRP = 1
    ISET = 1
  end if

  if (ISET == 1) then
    ! A new permutation was found, check and see if it was obtained previously
    do LPERM=1,NPERM
      if ((JATP(LPERM) == KATP) .and. (JASM(LPERM) == KASM) .and. (JBTP(LPERM) == KBTP) .and. (JBSM(LPERM) == KBSM)) exit
    end do
    if (LPERM > NPERM) then
      ! The permutation was new, add it to the list
      NPERM = NPERM+1
      JASM(NPERM) = KASM
      JBSM(NPERM) = KBSM
      JATP(NPERM) = KATP
      JBTP(NPERM) = KBTP
      if ((NPERM == 1) .or. ((NPERM >= 1) .and. (KSIGN == LSIGN))) then
        ISGN(NPERM) = 1
      else
        ISGN(NPERM) = -1
      end if
      LSIGN = KSIGN
      if ((NPERM == 1) .or. ((NPERM >= 1) .and. (KTRP == LTRP))) then
        ITRP(NPERM) = 0
      else
        ITRP(NPERM) = 1
      end if
      LTRP = KTRP
    end if
  end if
end do

! Should the block be trnasposed or scaled to return to initial form
ITRP(NPERM+1) = LTRP
ISGN(NPERM+1) = LSIGN
#ifdef _DEBUGPRINT_
write(u6,'(A,4I4)') ' Blocks obtained from IASM IBSM IATP IBTP ',IASM,IBSM,IATP,IBTP
write(u6,*)
write(u6,'(A)') ' JASM JBSM JATP JBTP Isgn Itrp'
write(u6,*)
do IPERM=1,NPERM
  write(u6,'(2x,6I4)') JASM(IPERM),JBSM(IPERM),JATP(IPERM),JBTP(IPERM),ISGN(IPERM),ITRP(IPERM)
end do
#endif

end subroutine PRMBLK
