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

subroutine SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END)
! Generate a set of sigma blocks,
! The NBLOCK specified in IBLOCK starting from IBOFF,
! be more specific.
!
! The blocks are delivered in HCB
!
! The blocks are scaled and reformed to combination order
! If LUCBLK > 0, the blocks of C corresponding to IBLOCK
! are stored on LUCBLK
!
! CONSPA,CONSPB  added October 1996
! ICBAT_RES, ICBAT_INI, IBBAT_END added august 1997
!
! If ICBAT_RES == 1 then it as assumed that only
! Cbatches ICBAT_INI to ICBAT_END are stored on  LUC

use stdalloc, only: mma_allocate, mma_deallocate
use GLBBAS, only: VEC3
use hidscr, only: ZSCR, ZOCSTR => OCSTR, REO, Z
use Local_Arrays, only: CLBT, CLEBT, CI1BT, CIBT, CBLTP, Allocate_Local_Arrays, Deallocate_Local_Arrays
use strbas, only: NSTSO
use lucia_data, only: NGAS, IPHGAS
! Definition of c and sigma spaces
use CandS, only: ICSM, ICSPC, ISSPC
use lucia_data, only: MXSOOB, MXNTTS, ISMOST
use lucia_data, only: IPRCIX, IPRDIA
use lucia_data, only: IADVICE, ICJKAIB, IH0INSPC, IH0SPC, IOCPTSPC, ISIMSYM, IUSE_PH, LCSBLK, MOCAA, MXINKA, NPTSPC
use lucia_data, only: IDC, PSSIGN
use lucia_data, only: MXNSTR, IBSPGPFTP, ISPGPFTP, MAX_STR_OC_BLK, MAX_STR_SPGP, MNHL, NELFGP, NELFSPGP, NHLFSPGP, NSTFSMSPGP
use lucia_data, only: IDISK
use lucia_data, only: NSMOB
use lucia_data, only: I12, IPART, IPERTOP, I_RES_AB
use lucia_data, only: MXTSOB, NTOOB, NOCOB, IOBPTS, ITSOB, NOBPTS
use lucia_data, only: NOCTYP
use lucia_data, only: NELEC
use lucia_data, only: MXPOBS, MXPNGAS, MXPNSMST
use csm_data, only: NSMST, NSMDX, NSMSX
use csm_data, only: ADSXA, SXDXSX
#ifdef _DEBUGPRINT_
use lucia_data, only: ICISTR
use Definitions, only: u6
#endif

implicit none
! =====
! Input
! =====
! Sigma blocks require
integer NBLOCK, IBOFF, LUC, IRESTRICT, LUCBLK, ICBAT_RES, ICBAT_INI, ICBAT_END
integer IBLOCK(8,*)

real*8 CB(*), HCB(*)
integer, allocatable :: CONSPA(:), CONSPB(:)
real*8, allocatable :: INSCR(:), INSCR2(:)
integer, allocatable :: STSTS(:), STSTD(:)
integer, allocatable, target :: CIOIO(:), SIOIO(:)
integer, pointer :: SCIOIO(:)
integer, allocatable :: I1(:), I2(:), I3(:), I4(:)
real*8, allocatable :: XI1S(:), XI2S(:), XI3S(:), XI4S(:)
real*8, allocatable :: LSCLFAC(:)
integer, allocatable :: SVST(:)
integer, allocatable :: H0SPC(:)
integer, external :: IMNMX
integer NTEST, IATP, IBTP, IATPM1, IBTPM1, IATPM2, IBTPM2, NOCTPA, NOCTPB, IOCTPA, IOCTPB, NAEL, NBEL, MXSTBL0, MAXA, MAXA0, &
        MAXA1, MAXB, MAXB0, MAXB1, MXSTBL, MAXI, MAXK, IOBTP, IOBSM, LSCR1, INTSCR, LSCR2, MAXIK, LSCR3, NTTS, LZSCR, LZ, K12, &
        I1234, IDOH2, MXADKBLK, MXADKBLK_AS, MXCIJA, MXCIJAB, MXCIJB, MXCJ, MXCJ_ALLSYM, MXSXBL, MXSXST, MX_NSPII

NTEST = 0
if (LUCBLK > 0) IDISK(LUCBLK) = 0

! Info for this internal space
! type of alpha and beta strings
IATP = 1
IBTP = 2
! alpha and beta strings with an electron removed
IATPM1 = 3
IBTPM1 = 4
! alpha and beta strings with two electrons removed
IATPM2 = 5
IBTPM2 = 6

! Number of supergroups
NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
! Offset for supergroups
IOCTPA = IBSPGPFTP(IATP)
IOCTPB = IBSPGPFTP(IBTP)

NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)

! connection matrices for supergroups

call mma_allocate(CONSPA,NOCTPA**2,Label='CONSPA')
call mma_allocate(CONSPB,NOCTPB**2,Label='CONSPB')
!    SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
call SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,CONSPA,IPRCIX)
call SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,CONSPB,IPRCIX)

! string sym, string sym => sx sym
! string sym, string sym => dx sym
call mma_allocate(STSTS,NSMST**2,Label='STSTS')
call mma_allocate(STSTD,NSMST**2,Label='STSTD')
call STSTSM(STSTS,STSTD,NSMST)
! Largest block of strings in zero order space
MXSTBL0 = MXNSTR
! Largest number of strings of given symmetry and type
MAXA = 0
MAXA0 = IMNMX(NSTSO(IATP)%I,NSMST*NOCTYP(IATP),2)
MAXA = max(MAXA,MAXA0)
if (NAEL >= 1) then
  MAXA1 = IMNMX(NSTSO(IATPM1)%I,NSMST*NOCTYP(IATPM1),2)
  MAXA = max(MAXA,MAXA1)
  if (NAEL >= 2) then
    MAXA1 = IMNMX(NSTSO(IATPM2)%I,NSMST*NOCTYP(IATPM2),2)
    MAXA = max(MAXA,MAXA1)
  end if
end if

MAXB = 0
MAXB0 = IMNMX(NSTSO(IBTP)%I,NSMST*NOCTYP(IBTP),2)
MAXB = max(MAXB,MAXB0)
if (NBEL >= 1) then
  MAXB1 = IMNMX(NSTSO(IBTPM1)%I,NSMST*NOCTYP(IBTPM1),2)
  MAXB = max(MAXB,MAXB1)
  if (NBEL >= 2) then
    MAXB1 = IMNMX(NSTSO(IBTPM2)%I,NSMST*NOCTYP(IBTPM2),2)
    MAXB = max(MAXB,MAXB1)
  end if
end if
MXSTBL = max(MAXA,MAXB)
#ifdef _DEBUGPRINT_
if (IPRCIX >= 3) write(u6,*) ' Largest block of strings with given symmetry and type',MXSTBL
#endif
! Largest number of resolution strings and spectator strings
! that can be treated simultaneously
MAXI = min(MXINKA,MXSTBL)
MAXK = min(MXINKA,MXSTBL)
! Largest active orbital block belonging to given type and symmetry
MXTSOB = 0
do IOBTP=1,NGAS
  do IOBSM=1,NSMOB
    MXTSOB = max(MXTSOB,NOBPTS(IOBTP,IOBSM))
  end do
end do
!write(u6,*) ' MXTSOB = ',MXTSOB
! Local scratch arrays for blocks of C and sigma
!if (ISIMSYM == 0) then
LSCR1 = MXSOOB
!else
!  LSCR1 = MXSOOB_AS
!end if
LSCR1 = max(LSCR1,LCSBLK)
#ifdef _DEBUGPRINT_
if (IPRCIX >= 3) write(u6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
#endif
! SCRATCH space for integrals
! A 4 index integral block with four indices belonging OS class
INTSCR = max(MXTSOB**4,NTOOB**2)
#ifdef _DEBUGPRINT_
if (IPRCIX >= 3) write(u6,*) ' Integral scratch space ',INTSCR
#endif
call mma_allocate(INSCR,INTSCR,Label='INSCR')
call mma_allocate(INSCR2,INTSCR,Label='INSCR2')
! Arrays giving allowed type combinations
call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')
call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
! Offsets for alpha and beta supergroups
IOCTPA = IBSPGPFTP(IATP)
IOCTPB = IBSPGPFTP(IBTP)
! sigma needed for MXRESC
call IAIBCM(ISSPC,SIOIO)
call IAIBCM(ICSPC,CIOIO)
! Arrays for additional symmetry operation
!if ((IDC == 3) .or. (IDC == 4)) then
!  call mma_allocate(SVST,NSMST,Label='SVST')
!  call SIGVST(SVST,NSMST)
!else
call mma_allocate(SVST,1,Label='SVST')
!end if

! scratch space for projected matrices and a CI block

! Scratch space for CJKAIB resolution matrices
! Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
if (ISSPC >= ICSPC) then
  SCIOIO => SIOIO
else
  SCIOIO => CIOIO
end if
call MXRESCPH(SCIOIO,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NSMST,NSTFSMSPGP,MXPNSMST,NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,NELFSPGP,MXCJ, &
              MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
#ifdef _DEBUGPRINT_
if (IPRCIX >= 3) then
  write(u6,*) 'SBLOCK : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM',MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM
  write(u6,*) 'SBLOCK : MXADKBLK ',MXADKBLK
  write(u6,*) ' MX_NSPII = ',MX_NSPII
end if
#endif
! For hardwired routines MXCIJAB is also used
LSCR2 = max(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
#ifdef _DEBUGPRINT_
if (IPRCIX >= 3) write(u6,*) ' Space for resolution matrices ',LSCR2

if (IPRCIX >= 3) write(u6,*) ' LSCR2 = ',LSCR2
#endif
! I assume memory was allocated for blocks, so
!
! vectors able to hold strings of given sym and type
MAXIK = max(MAXI,MAXK)
! I1 and Xi1s must also be able to hold largest st block
LSCR3 = max(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
call mma_allocate(I1,LSCR3,Label='I1')
call mma_allocate(I2,LSCR3,Label='I2')
call mma_allocate(I3,LSCR3,Label='I3')
call mma_allocate(I4,LSCR3,Label='I4')
call mma_allocate(XI1S,LSCR3,Label='XI1S')
call mma_allocate(XI2S,LSCR3,Label='XI2S')
call mma_allocate(XI3S,LSCR3,Label='XI3S')
call mma_allocate(XI4S,LSCR3,Label='XI4S')

! Some TTS arrays
NTTS = MXNTTS

! for partitioning of vector
call Allocate_Local_Arrays(NTTS,NSMST)
call ZBLTP(ISMOST(1,ICSM),NSMST,IDC,CBLTP,SVST)
! For scaling for each TTS block
call mma_allocate(LSCLFAC,8*NTTS,Label='LSCLFAC')

! Space for four blocks of string occupations and arrays of
! reordering arrays
! Also used to hold an NORB*NORB matrix
LZSCR = (max(NAEL,NBEL)+3)*(NOCOB+1)+2*NOCOB+NOCOB*NOCOB
LZ = (max(NAEL,NBEL)+2)*NOCOB
! Set up to two blocks for orbital conserving operator
K12 = 1
call mma_allocate(ZOCSTR,MAX_STR_OC_BLK,K12,Label='ZOCSTR')
I1234 = 2
call mma_allocate(REO,MAX_STR_SPGP,I1234,Label='REO')
call mma_allocate(Z,LZ,I1234,Label='Z')
call mma_allocate(ZSCR,LZSCR,Label='ZSCR')
! 4 arrays containing all strings of given sym. Dimension can  be
!   reduced to largest number of strings in alpha or beta.
!write(u6,*) ' SBLOCKS : MAX_STR_SPGP = ',MAX_STR_SPGP

if (I12 == 2) then
  IDOH2 = 1
else
  IDOH2 = 0
end if
! Place perturbation integrals over one body integrals
! INSERT_START
if (I12 == 2) then
  IDOH2 = 1
else
  IDOH2 = 0
end if

! Prepare for perturbation calculation

!if (IPERTOP /= 0) then
! Matrix specifying partiotioned spaces
call mma_allocate(H0SPC,NOCTPA*NOCTPB,Label='H0SPC')
call H0INTSPC(IH0SPC,NPTSPC,IOCPTSPC,NOCTPA,NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NGAS,MXPNGAS,H0SPC,NELFGP)
!  if (IH0SPC == 0) then
! Form of perturbation in subspace has not been defined,
! Use current IPART
IH0INSPC(1) = IPART
!  end if
!end if

! Jesper: Initializing ksvst
!KCJPA = 1 ! jwk-cleanup
!KSIPA = 1 ! jwk-cleanup
call SBLOCKS(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,VEC3,CIOIO,ISMOST(1,ICSM),CBLTP,NSTSO(IATP)%I,NSTSO(IBTP)%I,NAEL,IATP,NBEL,IBTP,IOCTPA, &
             IOCTPB,NOCTPA,NOCTPB,NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXK,MAXI,LSCR1,INSCR,VEC3,VEC3(1+LSCR2), &
             STSTS,STSTD,SXDXSX,ADSXA,NGAS,NELFSPGP,IDC,I1,XI1S,I2,XI2S,IDOH2,MXPOBS,SVST,PSSIGN,IPRDIA,LUC,ICJKAIB,VEC3, &
             VEC3(1+LSCR2),I3,XI3S,I4,XI4S,MXSXST,MXSXBL,MOCAA,CLBT,CLEBT,CI1BT,CIBT,IRESTRICT,CONSPA,CONSPB,LSCLFAC,IPERTOP, &
             IH0INSPC,H0SPC,ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,I_RES_AB,ISIMSYM,INSCR2)

! CALL SBLOCKS --> 91

if (IDC == 2) then
  ! reform
  call RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,1,NSMST,NSTSO(IATP)%I,NSTSO(IBTP)%I,IDC,PSSIGN,1,NTEST)
  ! scale
  call SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK,NSMST,NSTSO(IATP)%I,NSTSO(IBTP)%I,IDC,1,NTEST)
end if

if (LUCBLK > 0) call ITODS([-1],1,-1,LUCBLK)
! Eliminate local memory
call mma_deallocate(CONSPA)
call mma_deallocate(CONSPB)
call mma_deallocate(STSTS)
call mma_deallocate(STSTD)
call mma_deallocate(INSCR)
call mma_deallocate(INSCR2)
call mma_deallocate(CIOIO)
call mma_deallocate(SIOIO)
nullify(SCIOIO)
call Deallocate_Local_Arrays()
call mma_deallocate(I1)
call mma_deallocate(I2)
call mma_deallocate(I3)
call mma_deallocate(I4)
call mma_deallocate(XI1S)
call mma_deallocate(XI2S)
call mma_deallocate(XI3S)
call mma_deallocate(XI4S)
call mma_deallocate(LSCLFAC)
call mma_deallocate(ZOCSTR)
call mma_deallocate(REO)
call mma_deallocate(Z)
call mma_deallocate(ZSCR)
call mma_deallocate(SVST)
call mma_deallocate(H0SPC)

end subroutine SBLOCK
