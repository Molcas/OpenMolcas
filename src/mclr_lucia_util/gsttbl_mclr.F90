!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine GSTTBL_MCLR(C,CTT,IATP,IASM,IBTP,IBSM,NOCTPA,NOCTPB,NSASO,NSBSO,PSSIGN,ICOOSC,IDC,PLSIGN,LUC,SCR)
! obtain  determinant block (iatp iasm, ibtp ibsm)
! from vector packed in combination format according to IDC

use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: C(*), PSSIGN, PLSIGN
real(kind=wp), intent(_OUT_) :: CTT(*), SCR(*)
integer(kind=iwp), intent(in) :: IATP, IASM, IBTP, IBSM, NOCTPA, NOCTPB, NSASO(NOCTPA,*), NSBSO(NOCTPB,*), &
                                 ICOOSC(NOCTPA,NOCTPB,*), IDC, LUC
integer(kind=iwp) :: IBASE, IDUM(1), IMZERO, ISGVST, LBL, LCOMB, LDET, NAST, NBST, NCI, NCOL, NELMNT, NRI, NROW
real(kind=wp) :: PLSSGN, PSIGN

PSIGN = Zero ! dummy initialize

! =================
! Read in from disc
! =================
if (LUC /= 0) then
  call IFRMDS(IDUM,1,-1,LUC)
  LBL = IDUM(1)
  call FRMDSC_MCLR(SCR,LBL,-1,LUC,IMZERO)
  NAST = NSASO(IATP,IASM)
  NBST = NSBSO(IBTP,IBSM)
  ! ISGVST is undefined
  if (LBL /= 0) call SDCMRF_MCLR_2(CTT,SCR,IATP,IBTP,IASM,IBSM,NAST,NBST,IDC,PSSIGN,PLSIGN,ISGVST,LDET,LCOMB)
  ! ISGVST and PLSIGN missing to make it work for IDC = 3,4
else
  ! ===============
  ! Pack out from C
  ! ===============
  ! Permutation sign
  if (IDC == 2) then
    PSIGN = PSSIGN
  else if (IDC == 3) then
    PSIGN = PLSIGN
  else
    PSIGN = Zero
  end if
  PLSSGN = PLSIGN*PSSIGN
  ! check for different packing possibilities and unpack
  if ((IASM > IBSM) .or. (IDC == 1) .or. ((IDC == 3) .and. (IASM >= IBSM))) then
    !************
    ! IASM > IBSM
    !************
    if (IDC < 4) then
      ! Simple copy
      IBASE = ICOOSC(IATP,IBTP,IASM)
      NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
      CTT(1:NELMNT) = C(IBASE:IBASE+NELMNT-1)
    else if (IDC == 4) then
      ! MLMS packed
      if (IATP > IBTP) then
        IBASE = ICOOSC(IATP,IBTP,IASM)
        NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
        CTT(1:NELMNT) = C(IBASE:IBASE+NELMNT-1)
      else if (IATP == IBTP) then
        IBASE = ICOOSC(IATP,IATP,IASM)
        NAST = NSASO(IATP,IASM)
        call TRIPK2_2(CTT,C(IBASE),NAST,NAST,PLSIGN*PSSIGN)
      else if (IATP < IBTP) then
        IBASE = ICOOSC(IBTP,IATP,IASM)
        NROW = NSASO(IBTP,IASM)
        NCOL = NSBSO(IATP,IBSM)
        call TRNSPS(NROW,NCOL,C(IBASE),CTT)
        NELMNT = NROW*NCOL
        CTT(1:NELMNT) = PLSIGN*PSSIGN*CTT(1:NELMNT)
      end if
    end if
  else if (IASM == IBSM) then
    !************
    ! IASM = IBSM
    !************
    if ((IATP > IBTP) .or. (IDC == 3)) then
      ! simple copying
      IBASE = ICOOSC(IATP,IBTP,IASM)
      NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
      CTT(1:NELMNT) = C(IBASE:IBASE+NELMNT-1)
    else if (IATP == IBTP) then
      ! expand triangular packed matrix
      IBASE = ICOOSC(IATP,IBTP,IASM)
      NAST = NSASO(IATP,IASM)
      call TRIPK2_2(CTT,C(IBASE),NAST,NAST,PSSIGN)
    else if (IATP < IBTP) then
      ! transpose ibtp iasm iatp ibsm block
      IBASE = ICOOSC(IBTP,IATP,IASM)
      NRI = NSASO(IBTP,IASM)
      NCI = NSBSO(IATP,IASM)
      call TRNSPS(NRI,NCI,C(IBASE),CTT)
      if (PSSIGN == -One) CTT(1:NRI*NCI) = -CTT(1:NRI*NCI)
    end if
  else if (IASM < IBSM) then
    !************
    ! IASM < IBSM
    !************
    ! transpose ibtp ibsm iatp iasm block
    if (IDC < 4) then
      IBASE = ICOOSC(IBTP,IATP,IBSM)
      NRI = NSASO(IBTP,IBSM)
      NCI = NSBSO(IATP,IASM)
      if (IDC == 2) then
        call TRNSPS(NRI,NCI,C(IBASE),CTT)
      else if (IDC == 3) then
        CTT(1:NRI*NCI) = C(IBASE:IBASE+NRI*NCI-1)
      end if
      if (PSIGN == -One) CTT(1:NRI*NCI) = -CTT(NRI*NCI)
    else if (IDC == 4) then
      if (IBTP > IATP) then
        IBASE = ICOOSC(IBTP,IATP,IBSM)
        NRI = NSASO(IBTP,IBSM)
        NCI = NSBSO(IATP,IASM)
        call TRNSPS(NRI,NCI,C(IBASE),CTT)
        if (PSSIGN == -One) CTT(1:NRI*NCI) = -CTT(1:NRI*NCI)
      else if (IBTP == IATP) then
        IBASE = ICOOSC(IBTP,IATP,IBSM)
        NRI = NSASO(IATP,IBSM)
        NCI = NSBSO(IATP,IASM)
        call TRIPK2_2(CTT,C(IBASE),NRI,NCI,PLSSGN)
        if (PLSIGN == -One) CTT(1:NRI*NCI) = -CTT(1:NRI*NCI)
      else if (IBTP < IATP) then
        IBASE = ICOOSC(IATP,IBTP,IBSM)
        NELMNT = NSASO(IATP,IBSM)*NSBSO(IBTP,IASM)
        if (PLSIGN == -One) then
          CTT(1:NELMNT) = -C(IBASE:IBASE+NELMNT-1)
        else
          CTT(1:NELMNT) = C(IBASE:IBASE+NELMNT-1)
        end if
      end if
    end if
  end if
end if

end subroutine GSTTBL_MCLR
