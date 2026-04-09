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
!
! DETERMINE BASE ADDRESSES
! DFTP        : OPEN SHELL DETERMINANTS OF PROTO TYPE
! CFTP        : BRANCHING DIAGRAMS FOR PROTO TYPES
! DTOC        : CSF-DET TRANSFORMATION FOR PROTO TYPES
! CONF_OCC(I) : SPACE FOR STORING NCNSM CONFIGURATION EXPANSIONS

use Data_Structures, only: Alloc1DiArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: LOFFI = 8**6, MXPCSM = 8, MXPICI = 30, MXPIRR = 8, MXPNSMST = 8, MXPOBS = 8, MXPORB = 500, &
                                MXPPTSPC = 20, MXPR4T = 10, MXPSTT = 2500, MXPTSOB = 35
integer(kind=iwp), parameter :: MXPNGAS = MXPR4T+6

integer(kind=iwp) :: I12, I1234S, I12S, I2ELIMINATED_IN_GAS(MXPNGAS), I34S, I_AM_OUT(MXPSTT), I_ELIMINATE_GAS, I_RES_AB, IADVICE, &
                     IB_CONF_OCC(MXPORB+1), IB_CONF_REO(MXPORB+1), IB_SD_FOR_OPEN(MXPORB+1), IBGPSTR(MXPNGAS), IBSO(MXPOBS), &
                     IBSPGPFTP(MXPSTT), ICISTR, ICJKAIB, ICMBSPC(MXPSTT,MXPICI), IDC, IDIAG, IDISK(100), &
                     IELIMINATED_IN_GAS(MXPNGAS), IGSFGP(MXPSTT), IGSOCC(MXPNGAS,2), IGSOCCX(MXPNGAS,2,MXPICI), &
                     IH0INSPC(MXPPTSPC), IH0SPC, IH1FORM, INGRP_VAL, ini_h0, IOBPTS(6+MXPR4T,MXPOBS), IPART, IPHGAS(MXPNGAS), &
                     IPRCIX, IPRDEN, IREFSM, IREOST(MXPORB), IREOTS(MXPORB), IRESTR, ISIMSYM, ISMFSO(MXPORB), ISMFTO(MXPORB), &
                     ISPGPFTP(MXPNGAS,MXPSTT), ISTAC(MXPSTT,2), ITOOBS(MXPOBS), kvec3_length = 0, LCMBSPC(MXPICI), LCSBLK, LUC, &
                     LUDIA, LUHC, LUMOUT, LUSC1, LUSC2, LUSC3, LUSC34, LUSC35, LUSC36, LUSC37, LUSC38, LUSC39, LUSC40, &
                     MAX_STR_OC_BLK, MAX_STR_SPGP, MAXOP, Memory_Needed_Lucia = 0, MINMAX_SM_GP(2,MXPSTT), MINOP, MNGSOC(MXPNGAS), &
                     MNHL, MOCAA, MS2, MULTS, MXGSOC(MXPNGAS), MXINKA, MXNSTR, MXNTTS, MXSOOB, MXTSOB, N_2ELIMINATED_GAS, &
                     N_ELIMINATED_BATCHES, N_ELIMINATED_GAS, NACOB, NACOBS(MXPOBS), NACTEL, NBINT1, NBINT2, NCISPC, NCMBSPC, &
                     NCONF_ALL_SYM, NCONF_PER_OPEN(MXPORB+1,MXPCSM), NCONF_PER_SYM(MXPCSM), nconf_tot, nCSF_HEXS, &
                     NCSF_PER_SYM(MXPCSM), NDEOB, NELEC(MXPSTT), NELFGP(MXPSTT), NELFSPGP(MXPNGAS,MXPSTT), NELFTP(MXPSTT), &
                     NELIS(4), NGAS, NGPSTR(MXPNGAS), NGRP, NGSOBT(MXPNGAS), NGSSH(MXPIRR,MXPNGAS), NHLFSPGP(MXPSTT), NINOB, &
                     NINOBS(MXPOBS), NIRREP, NMXOCCLS, NOBPT(6+MXPR4T), NOBPTS(6+MXPR4T,MXPOBS), NOCOB, NOCSF, NOCTYP(MXPSTT), &
                     NOINT, NORB1, NORB2, NORB3, NPCMCNF(MXPORB+1), NPCSCNF(MXPORB+1), NPDTCNF(MXPORB+1), NPTSPC, NROOT, &
                     NSD_PER_SYM(MXPCSM), NSMOB, NSPGPFTP(MXPSTT), NSTFGP(MXPSTT), NSTFSMGP(MXPNSMST,MXPSTT), &
                     NSTFSMSPGP(MXPNSMST,MXPSTT), NSTRKS(4), NSTTP, NSTTYP, NTOOB, NTOOBS(MXPOBS), NTSPGP
real(kind=wp) :: ECORE, ECORE_HEX, ECORE_ORIG, PSSIGN, TDENSI(3), TSIGMA(6), XISPSM(MXPCSM,MXPICI)
logical(kind=iwp) :: Sigma_on_disk = .false.
character(len=6) :: ENVIRO
type(Alloc1DiArray_Type) :: CONF_OCC(8), CONF_REO(8), ISTSO(MXPSTT), NSTSO(MXPSTT), OCCSTR(MXPSTT), PGINT1(MXPOBS), &
                            PGINT1A(MXPOBS), STREO(MXPSTT), STSTM(MXPSTT,2), Zmat(MXPSTT)
type(Alloc1DiArray_Type), target :: SDREO_I(8)
type(Alloc1DiArray_Type), allocatable :: REO_PTDT(:), Z_PTDT(:)
integer(kind=iwp), allocatable :: CFTP(:), CI1BT(:), CIBT(:), CLBT(:), CLEBT(:), DFTP(:), IBCONF_ALL_SYM_FOR_OCCLS(:), IOCLS(:), &
                                  ISMDFGP(:), ISMSCR(:), ISTSGP(:), KINH1(:), KINH1_NOCCSYM(:), LSM1(:), LSM2(:), NACTSYM(:), &
                                  NSTSGP(:), OCSTR(:,:), PINT1(:), PINT2(:), REO(:,:), SPGPAN(:), SPGPCR(:), Z(:,:), ZSCR(:)
integer(kind=iwp), allocatable, target :: CBLTP(:)
integer(kind=iwp), pointer :: SDREO(:)
real(kind=wp), allocatable :: DStmp(:), Dtmp(:), DTOC(:), INT1(:), INT1O(:), PAtmp(:), Pscr(:), Ptmp(:), RF1(:), RF2(:), RHO1(:), &
                              SIGMA_VEC(:), SRHO1(:), VEC3(:)
real(kind=wp), allocatable, target :: CI_VEC(:)

public :: Allocate_Local_Arrays, CBLTP, CFTP, CI1BT, CI_VEC, CIBT, CLBT, CLEBT, CONF_OCC, CONF_REO, Deallocate_Local_Arrays, DFTP, &
          DStmp, Dtmp, DTOC, ECORE, ECORE_HEX, ECORE_ORIG, ENVIRO, I12, I1234S, I12S, I2ELIMINATED_IN_GAS, I34S, I_AM_OUT, &
          I_ELIMINATE_GAS, I_RES_AB, IADVICE, IB_CONF_OCC, IB_CONF_REO, IB_SD_FOR_OPEN, IBCONF_ALL_SYM_FOR_OCCLS, IBGPSTR, IBSO, &
          IBSPGPFTP, ICISTR, ICJKAIB, ICMBSPC, IDC, IDIAG, IDISK, IELIMINATED_IN_GAS, IGSFGP, IGSOCC, IGSOCCX, IH0INSPC, IH0SPC, &
          IH1FORM, INGRP_VAL, ini_h0, INT1, INT1O, IOBPTS, IOCLS, IPART, IPHGAS, IPRCIX, IPRDEN, IREFSM, IREOST, IREOTS, IRESTR, &
          ISIMSYM, ISMDFGP, ISMFSO, ISMFTO, ISMSCR, ISPGPFTP, ISTAC, ISTSGP, ISTSO, ITOOBS, KINH1, KINH1_NOCCSYM, kvec3_length, &
          LCMBSPC, LCSBLK, LOFFI, LSM1, LSM2, LUC, LUDIA, LUHC, LUMOUT, LUSC1, LUSC2, LUSC3, LUSC34, LUSC35, LUSC36, LUSC37, &
          LUSC38, LUSC39, LUSC40, MAX_STR_OC_BLK, MAX_STR_SPGP, MAXOP, Memory_Needed_Lucia, MINMAX_SM_GP, MINOP, MNGSOC, MNHL, &
          MOCAA, MS2, MULTS, MXGSOC, MXINKA, MXNSTR, MXNTTS, MXPCSM, MXPIRR, MXPNGAS, MXPNSMST, MXPOBS, MXPORB, MXPSTT, MXPTSOB, &
          MXSOOB, MXTSOB, N_2ELIMINATED_GAS, N_ELIMINATED_BATCHES, N_ELIMINATED_GAS, NACOB, NACOBS, NACTEL, NACTSYM, NBINT1, &
          NBINT2, NCISPC, NCMBSPC, NCONF_ALL_SYM, NCONF_PER_OPEN, NCONF_PER_SYM, NCONF_TOT, NCSF_HEXS, NCSF_PER_SYM, NDEOB, NELEC, &
          NELFGP, NELFSPGP, NELFTP, NELIS, NGAS, NGPSTR, NGRP, NGSOBT, NGSSH, NHLFSPGP, NINOB, NINOBS, NIRREP, NMXOCCLS, NOBPT, &
          NOBPTS, NOCOB, NOCSF, NOCTYP, NOINT, NORB1, NORB2, NORB3, NPCMCNF, NPCSCNF, NPDTCNF, NPTSPC, NROOT, NSD_PER_SYM, NSMOB, &
          NSPGPFTP, NSTFGP, NSTFSMGP, NSTFSMSPGP, NSTRKS, NSTSGP, NSTSO, NSTTP, NSTTYP, NTOOB, NTOOBS, NTSPGP, OCCSTR, OCSTR, &
          PAtmp, PGINT1, PGINT1A, PINT1, PINT2, Pscr, PSSIGN, Ptmp, REO, REO_PTDT, RF1, RF2, RHO1, SDREO, SDREO_I, Sigma_on_disk, &
          SIGMA_VEC, SPGPAN, SPGPCR, SRHO1, STREO, STSTM, TDENSI, TSIGMA, VEC3, XISPSM, Z, Z_PTDT, Zmat, ZSCR

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
