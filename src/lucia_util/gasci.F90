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

subroutine GASCI(ISM,ISPC)
! CI optimization in GAS space number ISPC for symmetry ISM
!
! Jeppe Olsen, Winter of 1995

use CandS, only: ICSM, ICSPC, ISSM, ISSPC
use lucia_data, only: Allocate_Local_Arrays, CBLTP, CI_VEC, CI1BT, CIBT, CLBT, CLEBT, Deallocate_Local_Arrays, ECORE, ECORE_ORIG, &
                      I12, I_ELIMINATE_GAS, I_RES_AB, IADVICE, IBSPGPFTP, ICISTR, IDC, IDIAG, IDISK, IGSOCC, IH1FORM, IPHGAS, &
                      IPRCIX, IREFSM, IRESTR, kvec3_length, LCSBLK, LUDIA, LUSC1, MNHL, MXINKA, MXNSTR, MXNTTS, MXPNGAS, MXPNSMST, &
                      MXSOOB, NCSF_PER_SYM, NELEC, NELFSPGP, NELFTP, NGAS, NHLFSPGP, NIRREP, NOBPT, NOBPTS, NOCSF, NOCTYP, NSMOB, &
                      NSTFSMSPGP, NSTSO, PSSIGN, VEC3, XISPSM
#ifdef _DEBUGPRINT_
use lucia_data, only: ICMBSPC, IGSOCCX, LCMBSPC
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ISM, ISPC
integer(kind=iwp) :: IATP, IATPM1, IATPM2, IBTP, IBTPM1, IBTPM2, IOCCLS_ARR(1), IOCTPA, IOCTPB, LBLOCK, LSCR12, LSCR2, MAXA, &
                     MAXA1, MAXB, MAXB1, MAXK, MX_NSPII, MXADKBLK, MXADKBLK_AS, MXCIJA, MXCIJAB, MXCIJB, MXCJ, MXCJ_ALLSYM, &
                     MXSTBL, MXSTBL0, MXSXBL, NAEL, NBATCH, NBEL, NBLOCK, NDET, NEL, NOCCLS, NOCTPA, NOCTPB, NTTS, NVAR, ZERO_ARR(1)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IGAS, II, JGASSPC, JJGASSPC
#endif
real(kind=wp) :: SHIFT
integer(kind=iwp), allocatable :: CIOIO(:), SVST(:)
integer(kind=iwp), external :: IMNMX

!MXACJ = 0
!MXACIJ = 0
!MXAADST = 0
! Normal integrals accessed
IH1FORM = 1
I_RES_AB = 0
! CI not CC
! Not just number conserving part
!IH_OCC_CONS_TEST = 0
!if (IH_OCC_CONS_TEST == 1) then
!  write(u6,*) ' IH_OCC_CONS set to one in GASCI'
!  IH_OCC_CONS = 1
!end if

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ====================================='
write(u6,*) ' Control has been transferred to GASCI'
write(u6,*) ' ====================================='
write(u6,*)
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
#endif

NDET = int(XISPSM(ISM,ISPC))
!NEL = NELCI(ISPC)
NEL = 0
#ifdef _DEBUGPRINT_
write(u6,*) ' Number of determinants/combinations  ',NDET
#endif
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
call OCCLS(1,NOCCLS,IOCCLS_ARR,NEL,NGAS,IGSOCC(:,1),IGSOCC(:,2),0,ZERO_ARR,NOBPT)

if (NOCSF == 1) then
  NVAR = NDET
else
  !*JESPER : Addition start
  !NVAR = NCSASM(ISM)
  NVAR = NCSF_PER_SYM(ISM)
  !*JESPER : Addition end
end if
#ifdef _DEBUGPRINT_
write(u6,*) '  NVAR in GASCI ',NVAR
#endif
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
call Allocate_Local_Arrays(NTTS,NIRREP)
! Additional info required to construct partitioning
call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')

call IAIBCM(ISPC,CIOIO)
call mma_allocate(SVST,1,Label='SVST')
call ZBLTP(ISM,NIRREP,IDC,CBLTP,SVST)
call mma_deallocate(SVST)

! Batches  of C vector
call PART_CIV2(IDC,NSTSO(IATP)%A,NSTSO(IBTP)%A,NOCTPA,NOCTPB,NIRREP,CIOIO,ISM,NBATCH,CLBT,CLEBT,CI1BT,CIBT,0)
! Number of BLOCKS
NBLOCK = CI1BT(NBATCH)+CLBT(NBATCH)-1

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
  MAXA1 = IMNMX(NSTSO(IATPM1)%A,NIRREP*NOCTYP(IATPM1),2)
  !write(u6,*) ' MAXA1 1',MAXA1
  MAXA = max(MAXA,MAXA1)
  if (NAEL >= 2) then
    MAXA1 = IMNMX(NSTSO(IATPM2)%A,NIRREP*NOCTYP(IATPM2),2)
    MAXA = max(MAXA,MAXA1)
  end if
end if
MAXB = 0
if (NBEL >= 1) then
  MAXB1 = IMNMX(NSTSO(IBTPM1)%A,NIRREP*NOCTYP(IBTPM1),2)
  MAXB = max(MAXB,MAXB1)
  if (NBEL >= 2) then
    MAXB1 = IMNMX(NSTSO(IBTPM2)%A,NIRREP*NOCTYP(IBTPM2),2)
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
call MXRESCPH(CIOIO,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NIRREP,NSTFSMSPGP,MXPNSMST,NSMOB,MXPNGAS,NGAS,NOBPTS,MAXK,NELFSPGP,MXCJ,MXCIJA, &
              MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
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
   ! Note that CI_VEC is used as a scratch array here and below.
  call GASDIAT(CI_VEC,LUDIA,SHIFT,ICISTR,I12,CBLTP,NBLOCK,CIBT)

  if ((NOCSF == 1) .and. (ICISTR == 1)) then
    IDISK(LUDIA) = 0
    call TODSC(CI_VEC,NVAR,-1,LUDIA)
  end if
  if (IPRCIX >= 2) write(u6,*) ' Diagonal constructed'
else
  write(u6,*) ' Diagonal not calculated'
end if

call Deallocate_Local_Arrays()
call mma_deallocate(CIOIO)
call mma_deallocate(VEC3)

end subroutine GASCI
