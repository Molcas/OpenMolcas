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

subroutine PSTTBL_MCLR(C,CTT,IATP,IASM,IBTP,IBSM,NOCTPA,NOCTPB,NSASO,NSBSO,PSIGN,ICOOSC,IAC,IDC,LUHC,SCR)
! add(IAC = 1) or copy (IAC =2) determinant block (iatp iasm, ibtp ibsm
! to vector packed in combination format
! iatp,iasm, ibtp,ibsm is assumed to be allowed combination block
!
! Combination type is defined by IDC
! IAC = 2  does not work for LUHC <= 0 !

use Constants, only: One

implicit real*8(A-H,O-Z)
dimension C(*), CTT(*), NSASO(NOCTPA,*), NSBSO(NOCTPB,*)
dimension ICOOSC(NOCTPA,NOCTPB,*)
dimension SCR(*)
dimension ISGVST(IBSM)
dimension IDUM(1)

! ======================
! Write directly to disc
! ======================

! Assumes complete block in,
! copies to lower half, scales  and write out.
if (LUHC /= 0) then
  NAST = NSASO(IATP,IASM)
  NBST = NSBSO(IBTP,IBSM)
  call SDCMRF_MCLR(CTT,SCR,1,IATP,IBTP,IASM,IBSM,NAST,NBST,IDC,PSIGN,PLSIGN,ISGVST,LDET,LCOMB)
  ! Note : PLSIGN and ISGVST missing in order to make it work for IDC=3,4
  IDUM(1) = LCOMB
  call ITODS(IDUM,1,-1,LUHC)
  call TODSC_MCLR(SCR,LCOMB,-1,LUHC)
else
  ! ==================
  ! Add to packed list
  ! ==================
  if ((IASM > IBSM) .or. (IDC == 1) .or. (IDC == 3)) then
    !*************
    !* IASM > IBSM
    !*************

    if (IDC < 4) then
      ! simple copying
      IBASE = ICOOSC(IATP,IBTP,IASM)
      NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
      if (IAC == 1) then
        call DaXpY_(NELMNT,One,CTT,1,C(IBASE),1)
      else if (IAC == 2) then
        call DCOPY_(NELMNT,CTT,1,C(IBASE),1)
      end if
    else if (IDC == 4) then
      if (IATP > IBTP) then
        IBASE = ICOOSC(IATP,IBTP,IASM)
        NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
        if (IAC == 1) then
          call DaXpY_(NELMNT,One,CTT,1,C(IBASE),1)
        else if (IAC == 2) then
          call DCOPY_(NELMNT,CTT,1,C(IBASE),1)
        end if
      else if (IATP == IBTP) then
        IBASE = ICOOSC(IATP,IBTP,IASM)
        NAST = NSASO(IATP,IASM)
        if (IAC == 1) then
          call PMPLFM(C(IBASE),CTT,NDIM)
        else
          call TRIPK2(CTT,C(IBASE),1,NAST,NAST,PSIGN)
        end if
      end if
    end if
  else if (IASM == IBSM) then
    !*************
    !* IASM = IBSM
    !*************
    if (IATP > IBTP) then
      ! simple copying
      IBASE = ICOOSC(IATP,IBTP,IASM)
      NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
      if (IAC == 1) then
        call DaXpY_(NELMNT,One,CTT,1,C(IBASE),1)
      else if (IAC == 2) then
        call DCOPY_(NELMNT,CTT,1,C(IBASE),1)
      end if
    else if (IATP == IBTP) then
      ! reform to triangular packed matrix
      IBASE = ICOOSC(IATP,IBTP,IASM)
      NAST = NSASO(IATP,IASM)
      if (IAC == 1) then
        call PMPLFM(C(IBASE),CTT,NAST)
      else
        call TRIPK2(CTT,C(IBASE),1,NAST,NAST,PSIGN)
      end if
    end if
  end if
end if

end subroutine PSTTBL_MCLR
