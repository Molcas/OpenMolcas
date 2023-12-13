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
! Copyright (C) 1989,2003, Jeppe Olsen                                 *
!***********************************************************************

subroutine CNHCN(ICNL,ITPL,ICNR,ITPR,CNHCNM,SCR,NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NORB,TUVX,IPREXH,ExFac,IREOTS)
! Obtain Hamiltonian matrix over CSF's of configurations ICNL,ICNR
!
! Jeppe Olsen, Summer of '89
!              IREOTS added August 2003

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICNL(*), ITPL, ICNR(*), ITPR, NAEL, NBEL, IPRODT(*), NORB, IREOTS(NORB)
real(kind=wp), intent(_OUT_) :: CNHCNM(*), SCR(*)
real(kind=wp), intent(in) :: ECORE, ONEBOD(NORB,NORB), DTOC(*), TUVX(*), ExFac
integer(kind=iwp), intent(inout) :: IPREXH
integer(kind=iwp) :: IPL, IPR, JCSF, JDET, KLCHD, KLDHD, KLDTLA, KLDTLB, KLDTRA, KLDTRB, KLFREE, KLISL, KLISR, NCSFL, NCSFR, &
                     NDETL, NDETR, NTEST
real(kind=wp) :: PSIGN
#include "spinfo.fh"
#include "ciinfo.fh"

call CNHCN_INTERNAL(SCR)

! This is to allow type punning without an explicit interface
contains

subroutine CNHCN_INTERNAL(SCR)

  real(kind=wp), target :: SCR(*)
  integer(kind=iwp), pointer :: iSCRla(:), iSCRlb(:), iSCRra(:), iSCRrb(:), iSCRf(:)
  integer(kind=iwp) :: JTYP

  NTEST = 0

  ! 1 : Obtain determinants corresponding to configurations

  NDETL = NDTFTP(ITPL)
  NDETR = NDTFTP(ITPR)

  NCSFL = NCSFTP(ITPL)
  NCSFR = NCSFTP(ITPR)

  KLFREE = 1

  KLDTLA = KLFREE
  KLFREE = KLFREE+NDETL*NAEL

  KLDTLB = KLFREE
  KLFREE = KLFREE+NDETL*NBEL

  KLISL = KLFREE
  KLFREE = KLFREE+NDETL

  KLDTRA = KLFREE
  KLFREE = KLFREE+NDETR*NAEL

  KLDTRB = KLFREE
  KLFREE = KLFREE+NDETR*NBEL

  KLISR = KLFREE
  KLFREE = KLFREE+NDETR

  call c_f_pointer(c_loc(SCR(KLDTLA)),iSCRla,[1])
  call c_f_pointer(c_loc(SCR(KLDTLB)),iSCRlb,[1])
  call c_f_pointer(c_loc(SCR(KLFREE)),iSCRf,[1])
  call CNFSTR(ICNL,ITPL,iSCRla,iSCRlb,NORB,NAEL,NBEL,NDETL,IPRODT,iSCRf,SCR(KLISL),IPREXH)
  nullify(iSCRla,iSCRlb,iSCRf)
  call c_f_pointer(c_loc(SCR(KLDTRA)),iSCRra,[1])
  call c_f_pointer(c_loc(SCR(KLDTRB)),iSCRrb,[1])
  call c_f_pointer(c_loc(SCR(KLFREE)),iSCRf,[1])
  call CNFSTR(ICNR,ITPR,iSCRra,iSCRrb,NORB,NAEL,NBEL,NDETR,IPRODT,iSCRf,SCR(KLISR),IPREXH)
  nullify(iSCRra,iSCRrb,iSCRf)

  ! Hamiltonian matrix over determinants of this configuration

  KLDHD = KLFREE
  KLFREE = KLFREE+NDETL*NDETR
  call COMBINATIONS(ICOMBI,PSIGN)
  call c_f_pointer(c_loc(SCR(KLDTLA)),iSCRla,[1])
  call c_f_pointer(c_loc(SCR(KLDTLB)),iSCRlb,[1])
  call c_f_pointer(c_loc(SCR(KLDTRA)),iSCRra,[1])
  call c_f_pointer(c_loc(SCR(KLDTRB)),iSCRrb,[1])
  call c_f_pointer(c_loc(SCR(KLFREE)),iSCRf,[1])
  call DIHDJ_MOLCAS(iSCRla,iSCRlb,NDETL,iSCRra,iSCRrb,NDETR,NAEL,NBEL,iSCRf,NORB,ONEBOD,SCR(KLDHD),0,0,ECORE,ICOMBI,PSIGN,0,TUVX, &
                    ExFac,IREOTS)
  nullify(iSCRla,iSCRlb,iSCRra,iSCRrb,iSCRf)

  ! Transform matrix to CSF basis

  ! : sign changes
  call DGMM2_MOLCAS(SCR(KLDHD),SCR(KLISL),1,NDETL,NDETR)
  call DGMM2_MOLCAS(SCR(KLDHD),SCR(KLISR),2,NDETL,NDETR)
  IPL = 1
  do JTYP=1,ITPL-1
    JCSF = NCSFTP(JTYP)
    JDET = NDTFTP(JTYP)
    IPL = IPL+JCSF*JDET
  end do
  if (ITPR == ITPL) then
    IPR = IPL
  else
    IPR = 1
    do JTYP=1,ITPR-1
      JCSF = NCSFTP(JTYP)
      JDET = NDTFTP(JTYP)
      IPR = IPR+JCSF*JDET
    end do
  end if

  KLCHD = KLFREE
  KLFREE = KLFREE+NCSFL*NDETR
  call MATML4(SCR(KLCHD),DTOC(IPL),SCR(KLDHD),NCSFL,NDETR,NDETL,NCSFL,NDETL,NDETR,1)
  call MATML4(CNHCNM,SCR(KLCHD),DTOC(IPR),NCSFL,NCSFR,NCSFL,NDETR,NDETR,NCSFR,0)

  if (NTEST >= 20) then
    write(u6,*) ' CSF-Hamiltonian matrix between two configurations'
    call WRTMAT(CNHCNM,NCSFL,NCSFR,NCSFL,NCSFR)
  end if

  return

end subroutine CNHCN_INTERNAL

end subroutine CNHCN
