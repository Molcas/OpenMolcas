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

subroutine CNHCN2(ICNL,ITPL,ICNR,ITPR,CNHCNM,SCR,NEL,NAEL,NBEL,INTSPC,NINOC,ECORE,IPRODT,DTOC,NORB,ICOMBI,PSSIGN,NTERMS,NDIF0, &
                  NDIF1,NDIF2)
! Obtain Hamilton matrix over CSFs of configurations ICNL,ICNR
!
! Jeppe Olsen, Summer of 89
!
! Modified for LUCIA, September 1993

implicit none
integer ITPL, ITPR, NEL, NAEL, NBEL, INTSPC, NINOC
real*8 ECORE
integer NORB, ICOMBI
real*8 PSSIGN
integer NTERMS, NDIF0, NDIF1, NDIF2
! Specific input
integer ICNL(*), ICNR(*)
! General input
integer IPRODT(*)
real*8 DTOC(*)
! Scratch
real*8 SCR(*)
! Output
real*8 CNHCNM(*)
! Interface to LUCIA common blocks in order to access strings
integer IDUMMY(1)

call CNHCN2_INTERNAL(SCR)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(NINOC)
  call Unused_real(ECORE)
  call Unused_integer(NTERMS)
end if

! This is to allow type punning without an explicit interface
contains

subroutine CNHCN2_INTERNAL(SCR)

  use iso_c_binding
  use MCLR_Data, only: IASTFI, IBSTFI
  use Str_Info, only: Str
  use MCLR_Data, only: MINOP, NCPCNT, NDPCNT

  implicit none
  integer IAGRP, IBGRP, IOPL, IOPR, ICLL, ICLR, NDETL, NDETR, NCSFL, NCSFR, KLFREE, KLDTLA, KLDTLB, KLISL, KLDTRA, KLDTRB, KLISR, &
          KLDHD, KLCHD, KLROU, LDIHDJ, LCNFST, NDIFF, ISYM, IPL, JTYP, JCSF, JDET, IPR, LWORK
  real*8 ECOREP
  real*8, target :: SCR(*)
  integer, pointer :: iSCR(:), iSCRa(:), iSCRb(:), iSCRar(:), iSCRbr(:), iSCRn(:), iSCRnn(:)

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

    isym = 0       ! eaw
    ecorep = 0.0d0 ! eaw
    call c_f_pointer(c_loc(SCR(KLROU+NEL)),iSCRn,[1])
    call c_f_pointer(c_loc(SCR(KLROU+2*NEL)),iSCRnn,[1])
    call DIHDJ2_MCLR(iSCRa,iSCRb,NDETL,iSCRar,iSCRbr,NDETR,NAEL,NBEL,iSCRnn,LWORK,NORB,SCR(KLDHD),ISYM,0,ECOREP,ICOMBI,PSSIGN, &
                     Str(IAGRP)%OCSTR,Str(IBGRP)%OCSTR,Str(IAGRP)%OCSTR,Str(IBGRP)%OCSTR,0,IDUMMY,IDUMMY,IDUMMY,IDUMMY,iSCR,iSCRn, &
                     NDIF0,NDIF1,NDIF2)
    nullify(iSCR,iSCRa,iSCRb,iSCRar,iSCRbr,iSCRn,iSCRnn)

    ! Transform matrix to CSF basis

    ! sign changes
    call DGMM2(SCR(KLDHD),SCR(KLDHD),SCR(KLISL),1,NDETL,NDETR)
    call DGMM2(SCR(KLDHD),SCR(KLDHD),SCR(KLISR),2,NDETL,NDETR)
    IPL = 1
    do JTYP=1,ITPL-1
      JCSF = NCPCNT(JTYP)
      JDET = NDPCNT(JTYP)
      IPL = IPL+JCSF*JDET
    end do
    if (ITPR == ITPL) then
      IPR = IPL
    else
      IPR = 1
      do JTYP=1,ITPR-1
        JCSF = NCPCNT(JTYP)
        JDET = NDPCNT(JTYP)
        IPR = IPR+JCSF*JDET
      end do
    end if

    call MATML4(SCR(KLCHD),DTOC(IPL),SCR(KLDHD),NCSFL,NDETR,NDETL,NCSFL,NDETL,NDETR,1)
    call MATML4(CNHCNM,SCR(KLCHD),DTOC(IPR),NCSFL,NCSFR,NCSFL,NDETR,NDETR,NCSFR,0)

  else if (NDIFF > 2) then
    CNHCNM(1:NCSFL*NCSFR) = 0.0d0
  end if

  return

end subroutine CNHCN2_INTERNAL

end subroutine CNHCN2
