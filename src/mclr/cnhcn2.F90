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
! Copyright (C) 1989,1993, Jeppe Olsen                                 *
!***********************************************************************

subroutine CNHCN2(ICNL,ITPL,ICNR,ITPR,CNHCNM,SCR,NEL,NAEL,NBEL,INTSPC,IPRODT,DTOC,NORB,ICOMBI,PSSIGN,NDIF0,NDIF1,NDIF2)
! Obtain Hamilton matrix over CSFs of configurations ICNL,ICNR
!
! Jeppe Olsen, Summer of 89
!
! Modified for LUCIA, September 1993

use iso_c_binding, only: c_f_pointer, c_loc
use MCLR_Data, only: IASTFI, IBSTFI, MINOP, NCPCNT, NDPCNT
use Str_Info, only: Str
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICNL(*), ITPL, ICNR(*), ITPR, NEL, NAEL, NBEL, INTSPC, IPRODT(*), NORB, ICOMBI
real(kind=wp), intent(_OUT_) :: CNHCNM(*), SCR(*)
real(kind=wp), intent(in) :: DTOC(*), PSSIGN
integer(kind=iwp), intent(out) :: NDIF0, NDIF1, NDIF2
integer(kind=iwp) :: IAGRP, IBGRP, ICLL, ICLR, IDUMMY(1), IOPL, IOPR, IPL, IPR, ISYM, KLCHD, KLDHD, KLDTLA, KLDTLB, KLDTRA, &
                     KLDTRB, KLFREE, KLISL, KLISR, KLROU, LCNFST, LDIHDJ, NCSFL, NCSFR, NDETL, NDETR, NDIFF
real(kind=wp) :: ECOREP

call CNHCN2_INTERNAL(SCR)

return

! This is to allow type punning without an explicit interface
contains

subroutine CNHCN2_INTERNAL(SCR)

  real(kind=wp), target, intent(_OUT_) :: SCR(*)
  integer(kind=iwp), pointer :: iSCR(:), iSCRa(:), iSCRar(:), iSCRb(:), iSCRbr(:), iSCRn(:), iSCRnn(:)

  ! Length of SCR : 6 * NDET + NDET**2 + NDET*NCSF +
  !                 MAX((NDET*NEL + 2*NEL),4*NORB + 2*NEL)
  ! where the last line arises from the local memory in CNFSTR and DIHDJ,
  ! respectively

  IAGRP = IASTFI(INTSPC)
  IBGRP = IBSTFI(INTSPC)

  ! 1: Obtain determinants corresponding to configurations

  IOPL = ITPL-1+MINOP
  IOPR = ITPR-1+MINOP

  ICLL = (NEL-IOPL)/2
  ICLR = (NEL-IOPR)/2

  NDETL = NDPCNT(ITPL)
  NDETR = NDPCNT(ITPR)

  NCSFL = NCPCNT(ITPL)
  NCSFR = NCPCNT(ITPR)

  KLFREE = 1
  ! 6* Ndet for holding string numbers and signs
  KLDTLA = KLFREE
  KLFREE = KLFREE+NDETL

  KLDTLB = KLFREE
  KLFREE = KLFREE+NDETL

  KLISL = KLFREE
  KLFREE = KLFREE+NDETL

  KLDTRA = KLFREE
  KLFREE = KLFREE+NDETR

  KLDTRB = KLFREE
  KLFREE = KLFREE+NDETR

  KLISR = KLFREE
  KLFREE = KLFREE+NDETR
  ! for transforming to CSFs
  KLDHD = KLFREE
  KLFREE = KLFREE+NDETL*NDETR

  KLCHD = KLFREE
  KLFREE = KLFREE+NCSFL*NDETR
  ! Scratch space for routines called (DIHDJ, CNFSTR)
  KLROU = KLFREE
  LDIHDJ = 4*NORB+2*NEL
  LCNFST = max(NDETL,NDETR)*NEL+2*NEL
  KLFREE = KLROU+max(LDIHDJ,LCNFST)

  ! Prescreen for CNFs that do not interact

  call c_f_pointer(c_loc(SCR(1)),iSCR,[1])
  call CMP2CN(ICNL,ICLL,IOPL,ICNR,ICLR,IOPR,iSCR,NORB,NDIFF)
  nullify(iSCR)

  if (NDIFF <= 2) then
    ! Strings for determinants of these configurations
    call c_f_pointer(c_loc(SCR(KLROU)),iSCR,[1])
    call c_f_pointer(c_loc(SCR(KLDTLA)),iSCRa,[1])
    call c_f_pointer(c_loc(SCR(KLDTLB)),iSCRb,[1])
    call CNFSTR_MCLR(ICNL,ITPL,iSCRa,iSCRb,NORB,NAEL,NBEL,NDETL,IPRODT,IAGRP,IBGRP,iSCR,SCR(KLISL))
    call c_f_pointer(c_loc(SCR(KLDTRA)),iSCRar,[1])
    call c_f_pointer(c_loc(SCR(KLDTRB)),iSCRbr,[1])
    call CNFSTR_MCLR(ICNR,ITPR,iSCRar,iSCRbr,NORB,NAEL,NBEL,NDETR,IPRODT,IAGRP,IBGRP,iSCR,SCR(KLISR))

    ! Hamiltonian matrix over determinants of the configurations

    isym = 0      ! eaw
    ecorep = Zero ! eaw
    call c_f_pointer(c_loc(SCR(KLROU+NEL)),iSCRn,[1])
    call c_f_pointer(c_loc(SCR(KLROU+2*NEL)),iSCRnn,[1])
    call DIHDJ2_MCLR(iSCRa,iSCRb,NDETL,iSCRar,iSCRbr,NDETR,NAEL,NBEL,iSCRnn,NORB,SCR(KLDHD),ISYM,ECOREP,ICOMBI,PSSIGN, &
                     Str(IAGRP)%OCSTR,Str(IBGRP)%OCSTR,Str(IAGRP)%OCSTR,Str(IBGRP)%OCSTR,0,IDUMMY,IDUMMY,IDUMMY,IDUMMY,iSCR,iSCRn, &
                     NDIF0,NDIF1,NDIF2)
    nullify(iSCR,iSCRa,iSCRb,iSCRar,iSCRbr,iSCRn,iSCRnn)

    ! Transform matrix to CSF basis

    ! sign changes
    call DGMM2(SCR(KLDHD),SCR(KLISL),1,NDETL,NDETR)
    call DGMM2(SCR(KLDHD),SCR(KLISR),2,NDETL,NDETR)
    IPL = 1+sum(NCPCNT(1:ITPL-1)*NDPCNT(1:ITPL-1))
    if (ITPR == ITPL) then
      IPR = IPL
    else
      IPR = 1+sum(NCPCNT(1:ITPR-1)*NDPCNT(1:ITPR-1))
    end if

    call MATML4(SCR(KLCHD),DTOC(IPL),SCR(KLDHD),NCSFL,NDETR,NDETL,NCSFL,NDETL,NDETR,1)
    call MATML4(CNHCNM,SCR(KLCHD),DTOC(IPR),NCSFL,NCSFR,NCSFL,NDETR,NDETR,NCSFR,0)

  else if (NDIFF > 2) then
    CNHCNM(1:NCSFL*NCSFR) = Zero
  end if

  return

end subroutine CNHCN2_INTERNAL

end subroutine CNHCN2
