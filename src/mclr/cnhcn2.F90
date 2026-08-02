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

use MCLR_Data, only: IASTFI, IBSTFI, MINOP, NCPCNT, NDPCNT
use Str_Info, only: Str
use Constants, only: Zero
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICNL(*), ITPL, ICNR(*), ITPR, NEL, NAEL, NBEL, INTSPC, IPRODT(*), NORB, ICOMBI
real(kind=wp), intent(_OUT_) :: CNHCNM(*), SCR(*)
real(kind=wp), intent(in) :: DTOC(*), PSSIGN
integer(kind=iwp), intent(out) :: NDIF0, NDIF1, NDIF2
integer(kind=iwp) :: IAGRP, IBGRP, ICLL, ICLR, IOPL, IOPR, IPL, IPR, ISYM, KLCHD, KLDHD, KLFREE, KLISL, KLISR, LCNFST, NCSFL, &
                     NCSFR, NDETL, NDETR, NDIFF
real(kind=wp) :: ECOREP
integer(kind=iwp), allocatable :: iSCR(:), iSCRl(:,:), iSCRn(:,:), iSCRr(:,:), jWork(:,:)

! Length of SCR : 2 * NDET + NDET**2 + NDET*NCSF +
!                 (NDET*NEL) + 2*NEL
! where the last line arises from the local memory in CNFSTR respectively

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
! 2* Ndet for holding string numbers and signs
KLISL = KLFREE
KLFREE = KLFREE+NDETL

KLISR = KLFREE
KLFREE = KLFREE+NDETR
! for transforming to CSFs
KLDHD = KLFREE
KLFREE = KLFREE+NDETL*NDETR

KLCHD = KLFREE
KLFREE = KLFREE+NCSFL*NDETR

! Prescreen for CNFs that do not interact

LCNFST = max(NDETL,NDETR)*NEL+2*NEL
call mma_allocate(iSCR,LCNFST,Label='iSCR')
call mma_allocate(jWork,NORB,4,Label='jWork')

call CMP2CN(ICNL,ICLL,IOPL,ICNR,ICLR,IOPR,jWork,NORB,NDIFF)

if (NDIFF <= 2) then
  call mma_allocate(iSCRl,NDETL,2,Label='iSCRl')
  call mma_allocate(iSCRr,NDETR,2,Label='iSCRr')
  call mma_allocate(iSCRn,max(NAEL,NBEL),2,Label='iSCRn')
  ! Strings for determinants of these configurations
  call CNFSTR_MCLR(ICNL,ITPL,iSCRl(:,1),iSCRl(:,2),NORB,NAEL,NBEL,NDETL,IPRODT,IAGRP,IBGRP,iSCR,SCR(KLISL))
  call CNFSTR_MCLR(ICNR,ITPR,iSCRr(:,1),iSCRr(:,2),NORB,NAEL,NBEL,NDETR,IPRODT,IAGRP,IBGRP,iSCR,SCR(KLISR))

  ! Hamiltonian matrix over determinants of the configurations

  isym = 0      ! eaw
  ecorep = Zero ! eaw
  call DIHDJ2_MCLR(iSCRl(:,1),iSCRl(:,2),NDETL,iSCRr(:,1),iSCRr(:,2),NDETR,NAEL,NBEL,jWork,NORB,SCR(KLDHD),ISYM,ECOREP,ICOMBI, &
                   PSSIGN,Str(IAGRP)%OCSTR,Str(IBGRP)%OCSTR,Str(IAGRP)%OCSTR,Str(IBGRP)%OCSTR,iSCRn(:,1),iSCRn(:,2),NDIF0,NDIF1, &
                   NDIF2)

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

  call mma_deallocate(iSCRl)
  call mma_deallocate(iSCRr)
  call mma_deallocate(iSCRn)

else if (NDIFF > 2) then
  CNHCNM(1:NCSFL*NCSFR) = Zero
end if

call mma_deallocate(iSCR)
call mma_deallocate(jWork)

end subroutine CNHCN2
