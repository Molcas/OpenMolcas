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

subroutine GSTTBL(C,CTT,IATP,IASM,IBTP,IBSM,NOCTPA,NOCTPB,NSASO,NSBSO,PSSIGN,ICOOSC,IDC,PLSIGN,LUC,SCR,NSMST,ISCALE,SCLFAC)
!***********************************************************************
! Variables status:
! C     = input CI vector
! CTT   = output CI vector in SD format
! obtain  determinant block (iatp iasm, ibtp ibsm)
! from vector packed in combination format according to IDC
!
! If ISCALE = 1, the routine scales and returns the block
! in determinant normalization, and SCLFAC = 1.0
!
! If ISCALE = 0, the routine does not perform any overall
! scaling, and a scale factor is returned in SCLFAC
!
! IF ISCALE = 0, zero blocks are not set explicitly to zero,
! instead zero is returned in SCLFAC
!
! ISCALE, SCLFAC added May 97

use lucia_data, only: IDISK
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: C(*)
real(kind=wp), intent(_OUT_) :: CTT(*), SCR(*)
integer(kind=iwp), intent(in) :: IATP, IASM, IBTP, IBSM, NOCTPA, NOCTPB, NSMST, NSASO(NSMST,*), NSBSO(NSMST,*), &
                                 ICOOSC(NOCTPA,NOCTPB,*), IDC, LUC, ISCALE
real(kind=wp), intent(in) :: PSSIGN, PLSIGN
real(kind=wp), intent(inout) :: SCLFAC
integer(kind=iwp) :: IAMPACK, IBASE, IDUMMY(1), IMZERO, ISGVST(1), LBL, LCOMB, LDET, NAST, NBST, NCI, NELMNT, NO_ZEROING, NRI
real(kind=wp) :: PSIGN

!write(u6,*) ' GSTTBL,  IATP,IASM,IBTP,IBSM,ISCALE'
!write(u6,*) IATP,IASM,IBTP,IBSM,ISCALE
! =================
! Read in from disc
! =================
if (LUC /= 0) then
  call IDAFILE(LUC,2,IDUMMY,1,IDISK(LUC))
  LBL = IDUMMY(1)
  call IDAFILE(LUC,2,IDUMMY,1,IDISK(LUC))
  !write(u6,*) ' LBL = ',LBL
  if (ISCALE == 1) then
    call FRMDSC(SCR,LBL,-1,LUC,IMZERO,IAMPACK)
  else
    NO_ZEROING = 1
    call FRMDSC2(SCR,LBL,-1,LUC,IMZERO,IAMPACK,NO_ZEROING)
  end if

  if ((IMZERO == 1) .and. (ISCALE == 0)) then
    SCLFAC = Zero
  else
    NAST = NSASO(IASM,IATP)
    NBST = NSBSO(IBSM,IBTP)
    if (LBL /= 0) then
      call SDCMRF(CTT,SCR,IATP,IBTP,IASM,IBSM,NAST,NBST,IDC,PSSIGN,PLSIGN,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
    else
      SCLFAC = Zero
    end if
  end if

  !write(u6,*) ' ISCALE and SCLFAC on return in GSTTBL',ISCALE,SCLFAC

else
  ! ISGVST and PLSIGN missing to make it work for IDC = 3,4
  ! =================
  ! Pack out from C
  ! =================
  if (ISCALE == 0) then
    write(u6,*) ' GSTTBL : LUC = 0 and ISCALE = 0'
    write(u6,*) ' I will scale as normal'
    SCLFAC = One
  end if
  ! Permutation sign
  ! To get rid of annoying compiler warning
  PSIGN = Zero
  if (IDC == 2) then
    PSIGN = PSSIGN
  else if (IDC == 3) then
    PSIGN = PLSIGN
  end if
  !PLSSGN = PLSIGN*PSSIGN
  ! check for different packing possibilities and unpack
  if ((IASM > IBSM) .or. (IDC == 1) .or. ((IDC == 3) .and. (IASM >= IBSM))) then
    !*************
    !* IASM > IBSM
    !*************
    if (IDC < 4) then
      ! Simple copy
      IBASE = ICOOSC(IATP,IBTP,IASM)
      NELMNT = NSASO(IASM,IATP)*NSBSO(IBSM,IBTP)
      CTT(1:NELMNT) = C(IBASE:IBASE+NELMNT-1)
      !write(u6,*) ' simple copy IBASE NELMNT ',IBASE,NELMNT
      !call WRTMAT(CTT,NSASO(IASM,IATP),NSBSO(IBSM,IBTP),NSASO(IASM,IATP),NSBSO(IBSM,IBTP))
    !idc else if (IDC == 4) then
    !      !MLMS packed
    !      if (IATP > IBTP) then
    !        IBASE = ICOOSC(IATP,IBTP,IASM)
    !        NELMNT = NSASO(IASM,IATP)*NSBSO(IBSM,IBTP)
    !        CTT(1:NELMNT) = C(IBASE:IBASE+NELMNT-1)
    !      else if (IATP == IBTP) then
    !        IBASE = ICOOSC(IATP,IATP,IASM)
    !        NAST = NSASO(IASM,IATP)
    !        call TRIPK32(CTT,C(IBASE),NAST,NAST,PLSIGN*PSSIGN)
    !      else if (IATP < IBTP) then
    !        IBASE = ICOOSC(IBTP,IATP,IASM)
    !        NROW = NSASO(IASM,IBTP)
    !        NCOL = NSBSO(IBSM,IATP)
    !        call TRPMT3(C(IBASE),NROW,NCOL,CTT)
    !        NELMNT = NROW*NCOL
    !        CTT(1:NELMNT) = PLSIGN*PSSIGN*CTT(1:NELMNT)
    !idc   end if
    end if
  else if (IASM == IBSM) then
    !*************
    !* IASM = IBSM
    !*************
    if ((IATP > IBTP) .or. (IDC == 3)) then
      ! simple copying
      IBASE = ICOOSC(IATP,IBTP,IASM)
      NELMNT = NSASO(IASM,IATP)*NSBSO(IBSM,IBTP)
      CTT(1:NELMNT) = C(IBASE:IBASE+NELMNT-1)
    else if (IATP == IBTP) then
      ! expand triangular packed matrix
      IBASE = ICOOSC(IATP,IBTP,IASM)
      NAST = NSASO(IASM,IATP)
      call TRIPK32(CTT,C(IBASE),NAST,NAST,PSSIGN)
    else if (IATP < IBTP) then
      ! transpose ibtp iasm iatp ibsm block
      IBASE = ICOOSC(IBTP,IATP,IASM)
      NRI = NSASO(IASM,IBTP)
      NCI = NSBSO(IASM,IATP)
      call TRPMT3(C(IBASE),NRI,NCI,CTT)
      if (PSSIGN == -One) CTT(1:NRI*NCI) = -CTT(1:NRI*NCI)
    end if
  else if (IASM < IBSM) then
    !*************
    !* IASM < IBSM
    !*************
    ! transpose ibtp ibsm iatp iasm block
    if (IDC < 4) then
      IBASE = ICOOSC(IBTP,IATP,IBSM)
      NRI = NSASO(IBSM,IBTP)
      NCI = NSBSO(IASM,IATP)
      if (IDC == 2) then
        call TRPMT3(C(IBASE),NRI,NCI,CTT)
      !else if (IDC == 3) then
      !  CTT(1:NRI*NCI) = C(IBASE:IBASE+NRI*NCI-1)
      end if
      if (PSIGN == -One) CTT(1:NRI*NCI) = -CTT(1:NRI*NCI)
    !idc else if (IDC  == 4) then
    !      if (IBTP > IATP) then
    !        IBASE = ICOOSC(IBTP,IATP,IBSM)
    !        NRI = NSASO(IBSM,IBTP)
    !        NCI = NSBSO(IASM,IATP)
    !        call TRPMT3(C(IBASE),NRI,NCI,CTT)
    !        if (PSSIGN == -One) CTT(1:NRI*NCI) = -CTT(1:NRI*NCI)
    !      else if (IBTP == IATP) then
    !        IBASE = ICOOSC(IBTP,IATP,IBSM)
    !        NRI   = NSASO(IBSM,IATP)
    !        NCI   = NSBSO(IASM,IATP)
    !        call TRIPK32(CTT,C(IBASE),NRI,NCI,PLSSGN)
    !        if (PLSIGN == -One) CTT(1:NRI*NCI) = -CTT(1:NRI*NCI)
    !      else if (IBTP < IATP) then
    !        IBASE = ICOOSC(IATP,IBTP,IBSM)
    !        NELMNT = NSASO(IBSM,IATP)*NSBSO(IASM,IBTP)
    !        CTT(1:NELMNT) = C(IBASE:IBASE+NELMNT-1)
    !        if (PLSIGN == -One) CTT(1:NELMNT) = -CTT(1:NELMNT)
    !idc   end if
    end if
  end if
end if

end subroutine GSTTBL
