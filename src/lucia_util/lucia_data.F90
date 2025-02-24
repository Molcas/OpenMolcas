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

module Lucia_Data

! Note : MXPNGAS = MXPR4T+6 !!
!        Required in order to handle GAS and RAS within /LUCINP/
! MXPPTSPC : Largest allowed division of space for perturbation operator
!
! CLBT  : Length of each Batch (in blocks)
! CLEBT : Length of each Batch (in elements)
! CI1BT : Length of each block
! CIBT  : Info on each block
! CBLTP : BLock type for each symmetry

use Data_Structures, only: Alloc1DiArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: LOFFI = 8**6, MXPCSM = 100, MXPICI = 30, MXPIRR = 20, MXPNSMST = 8, MXPOBS = 20, MXPORB = 500, &
                                MXPPTSPC = 20, MXPR4T = 10, MXPSTT = 2500, MXPTSOB = 35
integer(kind=iwp), parameter :: MXPNGAS = MXPR4T+6

integer(kind=iwp) :: I12, I1234S, I12S, I2ELIMINATED_IN_GAS(MXPNGAS), I34S, I_AM_OUT(MXPSTT), I_ELIMINATE_GAS, I_RES_AB, IADVICE, &
                     IB_CONF_OCC(MXPORB+1), IB_CONF_REO(MXPORB+1), IB_SD_FOR_OPEN(MXPORB+1), IBCONF_ALL_SYM_FOR_OCCLS(MXPCSM), &
                     IBGPSTR(MXPNGAS), IBSO(MXPOBS), IBSPGPFTP(MXPSTT), ICISTR, ICJKAIB, ICMBSPC(MXPSTT,MXPICI), IDC, IDIAG, &
                     IDISK(100), IELIMINATED_IN_GAS(MXPNGAS), IGSFGP(MXPSTT), IGSOCC(MXPNGAS,2), IGSOCCX(MXPNGAS,2,MXPICI), &
                     IH0INSPC(MXPPTSPC), IH0SPC, IH1FORM, INGRP_VAL, IOBPTS(6+MXPR4T,MXPOBS), IPART, IPHGAS(MXPNGAS), IPRCIX, &
                     IPRDEN, IREFSM, IREOST(MXPORB), IREOTS(MXPORB), IRESTR, ISIMSYM, ISMFSO(MXPORB), ISMFTO(MXPORB), &
                     ISMOST(MXPCSM,MXPCSM), ISPGPFTP(MXPNGAS,MXPSTT), ISTAC(MXPSTT,2), ITOOBS(MXPOBS), LCMBSPC(MXPICI), LCSBLK, &
                     LUC, LUDIA, LUHC, LUMOUT, LUSC1, LUSC2, LUSC3, LUSC34, LUSC35, LUSC36, LUSC37, LUSC38, LUSC39, LUSC40, &
                     MAX_STR_OC_BLK, MAX_STR_SPGP, MAXOP, MINMAX_SM_GP(2,MXPSTT), MINOP, MNGSOC(MXPNGAS), MNHL, MOCAA, MS2, MULTS, &
                     MXGSOC(MXPNGAS), MXINKA, MXNSTR, MXNTTS, MXSB, MXSOOB, MXTSOB, N_2ELIMINATED_GAS, N_ELIMINATED_BATCHES, &
                     N_ELIMINATED_GAS, NACOB, NACOBS(MXPOBS), NACTEL, NBINT1, NBINT2, NCISPC, NCMBSPC, NCONF_ALL_SYM, &
                     NCONF_PER_OPEN(MXPORB+1,MXPCSM), NCONF_PER_SYM(MXPCSM), nconf_tot, nCSF_HEXS, NCSF_PER_SYM(MXPCSM), NDEOB, &
                     NELEC(MXPSTT), NELFGP(MXPSTT), NELFSPGP(MXPNGAS,MXPSTT), NELFTP(MXPSTT), NELIS(4), NGAS, NGPSTR(MXPNGAS), &
                     NGRP, NGSOBT(MXPNGAS), NGSSH(MXPIRR,MXPNGAS), NHLFSPGP(MXPSTT), NINOB, NINOBS(MXPOBS), NIRREP, NMXOCCLS, &
                     NOBPT(6+MXPR4T), NOBPTS(6+MXPR4T,MXPOBS), NOCOB, NOCSF, NOCTYP(MXPSTT), NOINT, NORB1, NORB2, NORB3, &
                     NPCMCNF(MXPORB+1), NPCSCNF(MXPORB+1), NPDTCNF(MXPORB+1), NPTSPC, NROOT, NSD_PER_SYM(MXPCSM), NSMOB, &
                     NSPGPFTP(MXPSTT), NSTFGP(MXPSTT), NSTFSMGP(MXPNSMST,MXPSTT), NSTFSMSPGP(MXPNSMST,MXPSTT), NSTRKS(4), NSTTP, &
                     NSTTYP, NTOOB, NTOOBS(MXPOBS), NTSPGP
real(kind=wp) :: ECORE, ECORE_HEX, ECORE_ORIG, PSSIGN, XISPSM(MXPCSM,MXPICI)
character(len=6) :: ENVIRO
type(Alloc1DiArray_Type) :: ISTSO(MXPSTT), NSTSO(MXPSTT), OCCSTR(MXPSTT), STREO(MXPSTT), STSTM(MXPSTT,2), Zmat(MXPSTT)
integer(kind=iwp), allocatable :: CI1BT(:), CIBT(:), CLBT(:), CLEBT(:), IOCLS(:), ISMDFGP(:), ISMSCR(:), ISTSGP(:), NACTSYM(:), &
                                  NSTSGP(:), OCSTR(:,:), REO(:,:), SPGPAN(:), SPGPCR(:), Z(:,:), ZSCR(:)
integer(kind=iwp), allocatable, target :: CBLTP(:)

public :: Allocate_Local_Arrays, CBLTP, CI1BT, CIBT, CLBT, CLEBT, Deallocate_Local_Arrays, ECORE, ECORE_HEX, ECORE_ORIG, ENVIRO, &
          I12, I1234S, I12S, I2ELIMINATED_IN_GAS, I34S, I_AM_OUT, I_ELIMINATE_GAS, I_RES_AB, IADVICE, IB_CONF_OCC, IB_CONF_REO, &
          IB_SD_FOR_OPEN, IBCONF_ALL_SYM_FOR_OCCLS, IBGPSTR, IBSO, IBSPGPFTP, ICISTR, ICJKAIB, ICMBSPC, IDC, IDIAG, IDISK, &
          IELIMINATED_IN_GAS, IGSFGP, IGSOCC, IGSOCCX, IH0INSPC, IH0SPC, IH1FORM, INGRP_VAL, IOBPTS, IOCLS, IPART, IPHGAS, IPRCIX, &
          IPRDEN, IREFSM, IREOST, IREOTS, IRESTR, ISIMSYM, ISMDFGP, ISMFSO, ISMFTO, ISMOST, ISMSCR, ISPGPFTP, ISTAC, ISTSGP, &
          ISTSO, ITOOBS, LCMBSPC, LCSBLK, LOFFI, LUC, LUDIA, LUHC, LUMOUT, LUSC1, LUSC2, LUSC3, LUSC34, LUSC35, LUSC36, LUSC37, &
          LUSC38, LUSC39, LUSC40, MAX_STR_OC_BLK, MAX_STR_SPGP, MAXOP, MINMAX_SM_GP, MINOP, MNGSOC, MNHL, MOCAA, MS2, MULTS, &
          MXGSOC, MXINKA, MXNSTR, MXNTTS, MXPCSM, MXPIRR, MXPNGAS, MXPNSMST, MXPOBS, MXPORB, MXPSTT, MXPTSOB, MXSB, MXSOOB, &
          MXTSOB, N_2ELIMINATED_GAS, N_ELIMINATED_BATCHES, N_ELIMINATED_GAS, NACOB, NACOBS, NACTEL, NACTSYM, NBINT1, NBINT2, &
          NCISPC, NCMBSPC, NCONF_ALL_SYM, NCONF_PER_OPEN, NCONF_PER_SYM, NCONF_TOT, NCSF_HEXS, NCSF_PER_SYM, NDEOB, NELEC, NELFGP, &
          NELFSPGP, NELFTP, NELIS, NGAS, NGPSTR, NGRP, NGSOBT, NGSSH, NHLFSPGP, NINOB, NINOBS, NIRREP, NMXOCCLS, NOBPT, NOBPTS, &
          NOCOB, NOCSF, NOCTYP, NOINT, NORB1, NORB2, NORB3, NPCMCNF, NPCSCNF, NPDTCNF, NPTSPC, NROOT, NSD_PER_SYM, NSMOB, &
          NSPGPFTP, NSTFGP, NSTFSMGP, NSTFSMSPGP, NSTRKS, NSTSGP, NSTSO, NSTTP, NSTTYP, NTOOB, NTOOBS, NTSPGP, OCCSTR, OCSTR, &
          PSSIGN, REO, SPGPAN, SPGPCR, STREO, STSTM, XISPSM, Z, Zmat, ZSCR

contains

subroutine Allocate_Local_Arrays(MXNTTS,NSMST)

  integer(kind=iwp), intent(in) :: MXNTTS, NSMST

  if (allocated(CLBT)) call Abend()
  call mma_allocate(CLBT,MXNTTS,Label='CLBT',safe='*')
  call mma_allocate(CLEBT,MXNTTS,Label='CLEBT',safe='*')
  call mma_allocate(CI1BT,MXNTTS,Label='CI1BT',safe='*')
  call mma_allocate(CIBT,8*MXNTTS,Label='CIBT',safe='*')
  call mma_allocate(CBLTP,NSMST,Label='CBLTP',safe='*')

end subroutine Allocate_Local_Arrays

subroutine Deallocate_Local_Arrays()

  call mma_deallocate(CLBT,safe='*')
  call mma_deallocate(CLEBT,safe='*')
  call mma_deallocate(CI1BT,safe='*')
  call mma_deallocate(CIBT,safe='*')
  call mma_deallocate(CBLTP,safe='*')

end subroutine Deallocate_Local_Arrays

end module Lucia_Data
