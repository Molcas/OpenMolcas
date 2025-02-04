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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

subroutine GASCI(ISM,ISPC,IPRNT)
! CI optimization in GAS space number ISPC for symmetry ISM
!
! Jeppe Olsen, Winter of 1995

use stdalloc, only: mma_allocate, mma_deallocate
! Note that CI_VEC is used as a scratch array here and below.
use GLBBAS, only: VEC3, SCR => CI_VEC
use Local_Arrays, only: CLBT, CLEBT, CI1BT, CIBT, CBLTP, Allocate_Local_Arrays, Deallocate_Local_Arrays
use strbas, only: NSTSO
use rasscf_lucia, only: kvec3_length
! module for communicating with sigma
use CandS, only: ICSM, ISSM, ICSPC, ISSPC
use lucia_data, only: NCSF_PER_SYM
use lucia_data, only: ECORE_ORIG, ECORE
use lucia_data, only: IGSOCC, IPHGAS, NGAS
use lucia_data, only: MXSOOB, MXNTTS, IDUMMY, ISMOST, NELCI, XISPSM
use lucia_data, only: LUDIA, LUSC1
use lucia_data, only: IPRCIX
use lucia_data, only: NOCSF, IDIAG, IRESTR, ICISTR, IADVICE, ISIMSYM, LCSBLK, MXINKA
use lucia_data, only: IREFSM, PSSIGN, IDC
use lucia_data, only: I_ELIMINATE_GAS, MXNSTR, IBSPGPFTP, MNHL, NELFSPGP, NELFTP, NHLFSPGP, NSTFSMSPGP
use lucia_data, only: IH1FORM, IH2FORM
use lucia_data, only: IDISK
use lucia_data, only: NSMOB
use lucia_data, only: I_RES_AB, I12
use lucia_data, only: NOBPT, NOBPTS
use lucia_data, only: NOCTYP
use lucia_data, only: NELEC
use lucia_data, only: MXPNGAS, MXPNSMST
#ifdef _DEBUGPRINT_
use lucia_data, only: LCMBSPC, ICMBSPC, IGSOCCX
#endif
use csm_data, only: NSMST
use Constants, only: Zero, Two
use Definitions, only: u6

implicit none
integer ISM, ISPC, IPRNT
integer IOCCLS_ARR(1), ZERO_ARR(1)
integer, allocatable :: CIOIO(:)
integer, allocatable :: SVST(:)
integer NTEST, NDET, IATP, IBTP, NEL, NOCCLS, LBLOCK, NOCTPA, NOCTPB, NTTS, NBLOCK, MXSTBL0, IATPM1, IBTPM1, IATPM2, IBTPM2, NAEL, &
        NBEL, MAXA, MAXA1, MAXB, MAXB1, MXSTBL, MAXK, IOCTPA, IOCTPB, MXCIJA, MXCIJB, MXSXBL, MXADKBLK, MXADKBLK_AS, LSCR2, &
        LSCR12, MXCIJAB, NVAR, MXCJ_ALLSYM, MX_NSPII, NBATCH, MXCJ
integer, external :: IFRMR
integer, external :: IMNMX
real*8 SHIFT
#ifdef _DEBUGPRINT_
integer IGAS, II, JJGASSPC, JGASSPC
#endif
! Should all parameters be tranfered to Molcas?
!parameter (IALL = 0)

NTEST = 1
NTEST = max(NTEST,IPRNT)
!MXACJ = 0
!MXACIJ = 0
!MXAADST = 0
! Normal integrals accessed
IH1FORM = 1
I_RES_AB = 0
IH2FORM = 1
! CI not CC
! Not just number conserving part
!IH_OCC_CONS_TEST = 0
!if (IH_OCC_CONS_TEST == 1) then
!  write(u6,*) ' IH_OCC_CONS set to one in GASCI'
!  IH_OCC_CONS = 1
!end if

#ifdef _DEBUGPRINT_
if (NTEST >= 20) then
  write(u6,*)
  write(u6,*) ' ====================================='
  write(u6,*) ' Control has been transferred to GASCI'
  write(u6,*) ' ====================================='
  write(u6,*)
end if
if (NTEST >= 5) then
  write(u6,'(A)') '  A few pertinent data :'
  write(u6,*)
  write(u6,'(A,I2)') '  CI space         ',ISPC
  write(u6,*)
  write(u6,*) ' Number of GAS spaces included ',LCMBSPC(ISPC)
  write(u6,'(A,10I3)') '  GAS spaces included           ',(ICMBSPC(II,ISPC),II=1,LCMBSPC(ISPC))
  write(u6,*)
  write(u6,*) ' Occupation constraints :'
  write(u6,*) '========================='
  write(u6,*)
  write(u6,*)
  do JJGASSPC=1,LCMBSPC(ISPC)
    JGASSPC = ICMBSPC(JJGASSPC,ISPC)
    write(u6,*) ' Gas space  Min acc. occupation Max acc. occupation'
    write(u6,*) ' =================================================='
    do IGAS=1,NGAS
      write(u6,'(3X,I2,13X,I3,16X,I3)') IGAS,IGSOCCX(IGAS,1,JGASSPC),IGSOCCX(IGAS,2,JGASSPC)
    end do
  end do
end if
#endif

NDET = int(XISPSM(ISM,ISPC))
NEL = NELCI(ISPC)
if (NTEST >= 20) write(u6,*) ' Number of determinants/combinations  ',NDET
if (NDET == 0) then
  write(u6,*) ' The number of determinants/combinations is zero.'
  write(u6,*) ' I am sure that fascinating discussions about'
  write(u6,*) ' the energy of such a wave function exists,'
  write(u6,*) ' but I am just a dumb program, so I will stop'
  write(u6,*)
  write(u6,*) ' GASCI : Vanishing number of parameters'
  !stop ' GASCI : Vanishing number of parameters'
  call SYSABENDMSG('lucia_util/gasci','User error','')
end if
! Transfer to CANDS
ICSM = ISM
ISSM = ISM
ICSPC = ISPC
ISSPC = ISPC
! Complete operator
I12 = 2
! Class info
!if (ICLSSEL == 1) then
! Number of occupation classes
IATP = 1
IBTP = 2
NEL = NELFTP(IATP)+NELFTP(IBTP)
ZERO_ARR(1) = 0
call OCCLS(1,NOCCLS,IOCCLS_ARR,NEL,NGAS,IGSOCC(1,1),IGSOCC(1,2),0,ZERO_ARR,NOBPT)
! and then the occupation classes

if (NOCSF == 1) then
  NVAR = NDET
else
  !*JESPER : Addition start
  !NVAR = NCSASM(ISM)
  NVAR = NCSF_PER_SYM(ISM)
  !*JESPER : Addition end
end if
if (IPRNT >= 5) write(u6,*) '  NVAR in GASCI ',NVAR
! Allocate memory for diagonalization
!if (ISIMSYM == 0) then
LBLOCK = MXSOOB
!else
!  LBLOCK = MXSOOB_AS
!end if
LBLOCK = max(LBLOCK,LCSBLK)
! JESPER : Should reduce I/O
LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
if (PSSIGN /= Zero) LBLOCK = int(Two*XISPSM(IREFSM,1))

! Information about block structure- needed by new PICO2 routine.
! Memory for partitioning of C vector
IATP = 1
IBTP = 2
NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
NTTS = MXNTTS
call Allocate_Local_Arrays(NTTS,NSMST)
! Additional info required to construct partitioning
call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')

call IAIBCM(ISPC,CIOIO)
call mma_allocate(SVST,1,Label='SVST')
call ZBLTP(ISMOST(1,ISM),NSMST,IDC,CBLTP,SVST)
call mma_deallocate(SVST)

! Batches  of C vector
call PART_CIV2(IDC,CBLTP,NSTSO(IATP)%I,NSTSO(IBTP)%I,NOCTPA,NOCTPB,NSMST,LBLOCK,CIOIO,ISMOST(1,ISM),NBATCH,CLBT,CLEBT,CI1BT,CIBT, &
               0,ISIMSYM)
! Number of BLOCKS
NBLOCK = IFRMR(CI1BT,1,NBATCH)+IFRMR(CLBT,1,NBATCH)-1

! Enabling the calculation of excited states in a new way. Lasse
! This can be realized for GAS1 to GASN (for the GAS version)
! Here we find which type should be eliminated in the diagonal and
! in the sigma vector calculation. This to satisfy RASSI.
! Insert if statement
if (I_ELIMINATE_GAS >= 1) call I_AM_SO_EXCITED(NBATCH,CIBT,CLBT,CI1BT)
! End of story for Lasse

call mma_deallocate(CLBT)
call mma_deallocate(CLEBT)
! Length of each block
call EXTRROW(CIBT,8,8,NBLOCK,CI1BT)
call mma_deallocate(CI1BT)

! Class divisions of dets
! (Well, in principle I am against class division, but I
!  realize that it is a fact so ...)

! If PICO2/SBLOCK are used, three blocks are used in PICO2, so...

! Largest block of strings in zero order space
MXSTBL0 = MXNSTR
! type of alpha and beta strings
IATP = 1
IBTP = 2
! alpha and beta strings with an electron removed
IATPM1 = 3
IBTPM1 = 4
! alpha and beta strings with two electrons removed
IATPM2 = 5
IBTPM2 = 6

NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)
! Largest number of strings of given symmetry and type
MAXA = 0
if (NAEL >= 1) then
  MAXA1 = IMNMX(NSTSO(IATPM1)%I,NSMST*NOCTYP(IATPM1),2)
  !write(u6,*) ' MAXA1 1',MAXA1
  MAXA = max(MAXA,MAXA1)
  if (NAEL >= 2) then
    MAXA1 = IMNMX(NSTSO(IATPM2)%I,NSMST*NOCTYP(IATPM2),2)
    MAXA = max(MAXA,MAXA1)
  end if
end if
MAXB = 0
if (NBEL >= 1) then
  MAXB1 = IMNMX(NSTSO(IBTPM1)%I,NSMST*NOCTYP(IBTPM1),2)
  MAXB = max(MAXB,MAXB1)
  if (NBEL >= 2) then
    MAXB1 = IMNMX(NSTSO(IBTPM2)%I,NSMST*NOCTYP(IBTPM2),2)
    MAXB = max(MAXB,MAXB1)
  end if
end if
MXSTBL = max(MAXA,MAXB,MXSTBL0)
if (IPRCIX >= 2) write(u6,*) ' Largest block of strings with given symmetry and type',MXSTBL
! Largest number of resolution strings and spectator strings
! that can be treated simultaneously
MAXK = min(MXINKA,MXSTBL)
! scratch space for projected matrices and a CI block

! Scratch space for CJKAIB resolution matrices
! Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
IOCTPA = IBSPGPFTP(IATP)
IOCTPB = IBSPGPFTP(IBTP)
call MXRESCPH(CIOIO,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NSMST,NSTFSMSPGP,MXPNSMST,NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,NELFSPGP,MXCJ, &
              MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
if (IPRCIX >= 2) then
  write(u6,*) 'GASCI  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
  write(u6,*) ' MXADKBLK, MXADKBLK_AS',MXADKBLK,MXADKBLK_AS
end if
!if (ISIMSYM == 1) then
!  MXCJ = MAX(MXCJ_ALLSYM,MX_NSPII)
!  MXADKBLK_AS = MXADKBLK
!end if
! Using hardwired routines, MXCIJAB also used
LSCR2 = max(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
if (IPRCIX >= 2) write(u6,*) ' Space for two resolution matrices ',2*LSCR2
LSCR12 = max(LBLOCK,2*LSCR2)
call mma_allocate(VEC3,LSCR12,Label='VEC3')
KVEC3_LENGTH = max(LSCR12,2*LBLOCK,KVEC3_LENGTH)

! CI diagonal - if required

if (IDIAG == 2) LUDIA = LUSC1
if (.not. ((IDIAG == 2) .and. (IRESTR == 1))) then
  if (ICISTR >= 2) IDISK(LUDIA) = 0
  I12 = 2
  SHIFT = ECORE_ORIG-ECORE
  call GASDIAT(SCR,LUDIA,SHIFT,ICISTR,I12,CBLTP,NBLOCK,CIBT)

  if ((NOCSF == 1) .and. (ICISTR == 1)) then
    IDISK(LUDIA) = 0
    call TODSC(SCR,NVAR,-1,LUDIA)
  end if
  if (IPRCIX >= 2) write(u6,*) ' Diagonal constructed'
else
  write(u6,*) ' Diagonal not calculated'
end if

IDUMMY = 1
call Deallocate_Local_Arrays()
call mma_deallocate(CIOIO)
call mma_deallocate(VEC3)

end subroutine GASCI
