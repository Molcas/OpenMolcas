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

use CandS, only: ICSM, ICSPC, ISSPC
use lucia_data, only: Allocate_Local_Arrays, CBLTP, CI1BT, CIBT, CLBT, CLEBT, Deallocate_Local_Arrays, I12, I_RES_AB, IADVICE, &
                      IBSPGPFTP, ICJKAIB, IDC, IDISK, IH0INSPC, IH0SPC, IPART, IPHGAS, ISPGPFTP, MAX_STR_OC_BLK, MAX_STR_SPGP, &
                      MNHL, MOCAA, MXINKA, MXNSTR, MXNTTS, MXPNGAS, MXPNSMST, MXTSOB, NELEC, NELFGP, NELFSPGP, NGAS, NHLFSPGP, &
                      NIRREP, NOBPTS, NOCOB, NOCTYP, NPTSPC, NSMOB, NSTFSMSPGP, NSTSO, NTOOB, OCSTR, PSSIGN, REO, VEC3, Z, ZSCR
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use lucia_data, only: ICISTR, IPRCIX
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NBLOCK, IBLOCK(8,*), IBOFF, LUC, IRESTRICT, LUCBLK, ICBAT_RES, ICBAT_INI, ICBAT_END
real(kind=wp), intent(inout) :: CB(*), HCB(*)
integer(kind=iwp) :: DUM(1), I1234, IATP, IATPM1, IATPM2, IBTP, IBTPM1, IBTPM2, IDOH2, INTSCR, IOCTPA, IOCTPB, K12, LSCR2, LSCR3, &
                     LZ, LZSCR, MAXA, MAXA0, MAXA1, MAXB, MAXB0, MAXB1, MAXI, MAXIK, MAXK, MX_NSPII, MXADKBLK, MXADKBLK_AS, &
                     MXCIJA, MXCIJAB, MXCIJB, MXCJ, MXCJ_ALLSYM, MXSTBL, MXSTBL0, MXSXBL, NAEL, NBEL, NOCTPA, NOCTPB, NTTS
integer(kind=iwp), allocatable :: CONSPA(:), CONSPB(:), H0SPC(:), I1(:), I2(:), I3(:), I4(:), SVST(:)
integer(kind=iwp), allocatable, target :: CIOIO(:), SIOIO(:)
integer(kind=iwp), pointer :: SCIOIO(:)
real(kind=wp), allocatable :: INSCR(:), LSCLFAC(:), XI1S(:), XI2S(:), XI3S(:), XI4S(:)
integer(kind=iwp), external :: IMNMX

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
!    SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON)
call SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,CONSPA)
call SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,CONSPB)

! Largest block of strings in zero order space
MXSTBL0 = MXNSTR
! Largest number of strings of given symmetry and type
MAXA = 0
MAXA0 = IMNMX(NSTSO(IATP)%A,NIRREP*NOCTYP(IATP),2)
MAXA = max(MAXA,MAXA0)
if (NAEL >= 1) then
  MAXA1 = IMNMX(NSTSO(IATPM1)%A,NIRREP*NOCTYP(IATPM1),2)
  MAXA = max(MAXA,MAXA1)
  if (NAEL >= 2) then
    MAXA1 = IMNMX(NSTSO(IATPM2)%A,NIRREP*NOCTYP(IATPM2),2)
    MAXA = max(MAXA,MAXA1)
  end if
end if

MAXB = 0
MAXB0 = IMNMX(NSTSO(IBTP)%A,NIRREP*NOCTYP(IBTP),2)
MAXB = max(MAXB,MAXB0)
if (NBEL >= 1) then
  MAXB1 = IMNMX(NSTSO(IBTPM1)%A,NIRREP*NOCTYP(IBTPM1),2)
  MAXB = max(MAXB,MAXB1)
  if (NBEL >= 2) then
    MAXB1 = IMNMX(NSTSO(IBTPM2)%A,NIRREP*NOCTYP(IBTPM2),2)
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
MXTSOB = max(0,maxval(NOBPTS(1:NGAS,1:NSMOB)))
!write(u6,*) ' MXTSOB = ',MXTSOB
! Local scratch arrays for blocks of C and sigma
!if (ISIMSYM == 0) then
!  LSCR1 = MXSOOB
!else
!  LSCR1 = MXSOOB_AS
!end if
!LSCR1 = max(LSCR1,LCSBLK)
#ifdef _DEBUGPRINT_
if (IPRCIX >= 3) write(u6,*) ' ICISTR ',ICISTR
#endif
! SCRATCH space for integrals
! A 4 index integral block with four indices belonging OS class
INTSCR = max(MXTSOB**4,NTOOB**2)
#ifdef _DEBUGPRINT_
if (IPRCIX >= 3) write(u6,*) ' Integral scratch space ',INTSCR
#endif
call mma_allocate(INSCR,INTSCR,Label='INSCR')
! Arrays giving allowed type combinations
call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')
call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
! Offsets for alpha and beta supergroups
IOCTPA = IBSPGPFTP(IATP)
IOCTPB = IBSPGPFTP(IBTP)
! sigma needed for MXRESCPH
call IAIBCM(ISSPC,SIOIO)
call IAIBCM(ICSPC,CIOIO)
! Arrays for additional symmetry operation
!if ((IDC == 3) .or. (IDC == 4)) then
!  call mma_allocate(SVST,NIRREP,Label='SVST')
!  call SIGVST(SVST,NIRREP)
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
call MXRESCPH(SCIOIO,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NIRREP,NSTFSMSPGP,MXPNSMST,NSMOB,MXPNGAS,NGAS,NOBPTS,MAXK,NELFSPGP,MXCJ,MXCIJA, &
              MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
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
call Allocate_Local_Arrays(NTTS,NIRREP)
call ZBLTP(ICSM,NIRREP,IDC,CBLTP,SVST)
! For scaling for each TTS block
call mma_allocate(LSCLFAC,8*NTTS,Label='LSCLFAC')

! Space for four blocks of string occupations and arrays of
! reordering arrays
! Also used to hold an NORB*NORB matrix
LZSCR = (max(NAEL,NBEL)+3)*(NOCOB+1)+2*NOCOB+NOCOB*NOCOB
LZ = (max(NAEL,NBEL)+2)*NOCOB
! Set up to two blocks for orbital conserving operator
K12 = 1
call mma_allocate(OCSTR,MAX_STR_OC_BLK,K12,Label='OCSTR')
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
call H0INTSPC(IH0SPC,NPTSPC,NOCTPA,NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NGAS,MXPNGAS,H0SPC,NELFGP)
!  if (IH0SPC == 0) then
! Form of perturbation in subspace has not been defined,
! Use current IPART
IH0INSPC(1) = IPART
!  end if
!end if

! Jesper: Initializing ksvst
!KCJPA = 1 ! jwk-cleanup
!KSIPA = 1 ! jwk-cleanup
call SBLOCKS(NBLOCK,IBLOCK(:,IBOFF),CB,HCB,VEC3,CIOIO,ICSM,NSTSO(IATP)%A,NSTSO(IBTP)%A,NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA, &
             NOCTPB,NIRREP,NSMOB,NOBPTS,MXPNGAS,MAXK,MAXI,INSCR,VEC3,VEC3(1+LSCR2),NGAS,NELFSPGP,IDC,I1,XI1S,I2,XI2S,IDOH2,SVST, &
             PSSIGN,LUC,ICJKAIB,VEC3,VEC3(1+LSCR2),I3,XI3S,I4,XI4S,MOCAA,CLBT,CLEBT,CI1BT,CIBT,IRESTRICT,CONSPA,CONSPB,LSCLFAC, &
             H0SPC,ICBAT_RES,ICBAT_INI,ICBAT_END,IPHGAS,I_RES_AB)

! CALL SBLOCKS --> 91

if (IDC == 2) then
  ! reform
  call RFTTS(HCB,CB,IBLOCK(:,IBOFF),NBLOCK,NIRREP,NSTSO(IATP)%A,NSTSO(IBTP)%A,IDC)
  ! scale
  call SCDTTS(HCB,IBLOCK(:,IBOFF),NBLOCK,NIRREP,NSTSO(IATP)%A,NSTSO(IBTP)%A,IDC)
end if

if (LUCBLK > 0) then
  DUM(1) = -1
  call ITODS(DUM,1,-1,LUCBLK)
end if
! Eliminate local memory
call mma_deallocate(CONSPA)
call mma_deallocate(CONSPB)
call mma_deallocate(INSCR)
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
call mma_deallocate(OCSTR)
call mma_deallocate(REO)
call mma_deallocate(Z)
call mma_deallocate(ZSCR)
call mma_deallocate(SVST)
call mma_deallocate(H0SPC)

end subroutine SBLOCK
