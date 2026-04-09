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

!#define _DEBUGPRINT_
subroutine SigmaVec(C,HC,kic)
! Outer routine for sigma vector generation
! RAS space
! This is a driver routine change in some rows just so it s
! should work in my mclr program. Of course I have to change
! the name on it but if you look in Jeppe's mv7old(relaci) you
! will find a subroutine that looks very much like this
! This routine is just a setup routine for memory etc

use Str_Info, only: CNSM, DTOC, ITYP_Dummy, nElec, NOCTYP, STR
use MCLR_Data, only: i12, IASTFI, IBSTFI, IBTSOB, ICISTR, IDC, IDIAG, ist, ITSOB, MAXI, MAXK, MNR1IC, MXINKA, MXR3IC, MXSB, &
                     MXSOOB, NOCSF, NOPART, NTSOB, PSSIGN, XISPSM
use CandS, only: ICSM, ICSPC, ISSM, ISSPC
use input_mclr, only: nIrrep, TimeDep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: C(*)
real(kind=wp), intent(_OUT_) :: HC(*)
integer(kind=iwp), intent(in) :: kic(2)
integer(kind=iwp) :: IATP, IATP1, IATP2, IBTP, IBTP1, IBTP2, IDOH2, idummy(1), IICOPY, IMNMX, INTSCR, JATP, JBTP, LLUC, LLUHC, &
                     LSCR1, LSCR12, LSCR2, LSCR3, LUC, LUHC, MAXA, MAXA0, MAXA1, MAXB, MAXB0, MAXB1, MAXE, MAXEL3, MAXIK, MAXPK, &
                     MOCAA, MXCIJA, MXCIJAB, MXCIJB, MXCJ, MXIJST, MXIJSTF, MXSTBL, MXSTBL0, MXSXBL, MXSXST, MXTSOB, NAEL, NBEL, &
                     NOCTPA, NOCTPB, NOOS, NSDET
integer(kind=iwp), allocatable :: CBLTP(:), CIOIO(:), I1(:), I2(:), I3(:), I4(:), OOS(:,:), SBLTP(:), SIOIO(:), SVST(:)
real(kind=wp), allocatable :: CB(:), INSCR(:), SB(:), XI1S(:), XI2S(:), XI3S(:), XI4S(:)
real(kind=wp), allocatable, target :: C2(:), CJRES(:), SIRES(:)
real(kind=wp), pointer :: pC2(:), pCJRES(:), pSIRES(:)

LUC = 0
LUHC = 0
IDUMMY(1) = 0
NSDET = nint(XISPSM(ISSM,ISSPC))

! The story of MV7 : I started out from nothing, absolutely zero,

HC(1:NSDET) = Zero

! Info for this internal space

IATP = IASTFI(ISSPC)
IBTP = IBSTFI(ISSPC)
JATP = IASTFI(ICSPC)
JBTP = IBSTFI(ICSPC)
if ((IATP /= JATP) .or. (IBTP /= JBTP)) then
  write(u6,*) ' My world is falling apart'
  write(u6,*) ' C and sigma belongs to different types of strings'
  write(u6,*) ' IATP IBTP JATP JBTP ',IATP,IBTP,JATP,JBTP
  call Abend()
end if
NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)

! Largest block of strings in zero order space
MAXA0 = IMNMX(Str(IATP)%NSTSO,nIrrep*NOCTYP(IATP),2)
MAXB0 = IMNMX(Str(IBTP)%NSTSO,nIrrep*NOCTYP(IBTP),2)
MXSTBL0 = max(MAXA0,MAXB0)
! Largest number of strings of given symmetry and type
MAXA = 0
if (NAEL >= 1) then
  MAXA1 = IMNMX(Str(IATP+1)%NSTSO,nIrrep*NOCTYP(IATP+1),2)
  MAXA = max(MAXA,MAXA1)
end if
if (NAEL >= 2) then
  MAXA1 = IMNMX(Str(IATP+2)%NSTSO,nIrrep*NOCTYP(IATP+2),2)
  MAXA = max(MAXA,MAXA1)
end if
MAXB = 0
if (NBEL >= 1) then
  MAXB1 = IMNMX(Str(IBTP+1)%NSTSO,nIrrep*NOCTYP(IBTP+1),2)
  MAXB = max(MAXB,MAXB1)
end if
if (NBEL >= 2) then
  MAXB1 = IMNMX(Str(IBTP+2)%NSTSO,nIrrep*NOCTYP(IBTP+2),2)
  MAXB = max(MAXB,MAXB1)
end if
MXSTBL = max(MAXA,MAXB)
#ifdef _DEBUGPRINT_
write(u6,*) ' Largest block of strings with given symmetry and type',MXSTBL
#endif
MAXI = min(MXINKA,MXSTBL)
MAXK = min(MXINKA,MXSTBL)
! Largest active orbital block belonging to given type and symmetry
MXTSOB = IMNMX(NTSOB,3*nIrrep,2)
! Local scratch arrays for blocks of C and sigma
LSCR1 = 0
if (ICISTR <= 2) then
  LSCR1 = MXSB
else if (ICISTR == 3) then
  LSCR1 = MXSOOB
end if
if (ICISTR == 1) then
  call mma_allocate(CB,LSCR1,Label='CB')
  call mma_allocate(SB,LSCR1,Label='SB')
end if
! SCRATCH space for integrals
! A 4 index integral block with four indices belonging OS class
INTSCR = MXTSOB**4

call mma_allocate(INSCR,INTSCR,Label='INSCR')
! Arrays giving allowed type combinations

call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')

! SIOIO,CIOIO [NOCTPA x NOCTPB]
! Allowed combinations of alpha and beta strings
! Sigma and CI vector respectively

! (RAS1 & RAS3 constrains)

call IAIBCM_MCLR(MNR1IC(ISSPC),MXR3IC(ISSPC),NOCTPA,NOCTPB,Str(IATP)%EL1,Str(IATP)%EL3,Str(IBTP)%EL1,Str(IBTP)%EL3,SIOIO)

call IAIBCM_MCLR(MNR1IC(ICSPC),MXR3IC(ICSPC),NOCTPA,NOCTPB,Str(IATP)%EL1,Str(IATP)%EL3,Str(IBTP)%EL1,Str(IBTP)%EL3,CIOIO)

! Arrays giving block type
call mma_allocate(SBLTP,nIrrep,Label='SBLTP')
call mma_allocate(CBLTP,nIrrep,Label='CBLTP')
! Arrays for additional symmetry operation

call mma_allocate(SVST,nIrrep,Label='SVST')

! scratch space for projected matrices and a CI block

! Scratch space for CJKAIB resolution matrices
! Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
if (NOPART == 0) then
  MAXPK = MAXK
else
  MAXPK = 0
end if

! Get memory requirements

! MXCJ:MXCIJA:MXCIJB:MXCIJAB:MXSXBL:MXIJST:MXIJSTF

IATP1 = min(IATP+1,ITYP_DUMMY)
IBTP1 = min(IbTP+1,ITYP_DUMMY)
IATP2 = min(IATP+2,ITYP_DUMMY)
IBTP2 = min(IbTP+2,ITYP_DUMMY)
call MXRESC(CIOIO,IATP,IBTP,NOCTPA,NOCTPB,nIrrep,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,Str(IATP1)%NSTSO,NOCTYP(IATP1),Str(IBTP1)%NSTSO, &
            NOCTYP(IBTP1),NTSOB,MAXpK,Str(IATP2)%NSTSO,NOCTYP(IATP2),Str(IBTP2)%NSTSO,NOCTYP(IBTP2),Str(IATP)%EL123, &
            Str(IBTP)%EL123,MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXIJST,MXIJSTF)

!xvectors able to hold strings of given sym and type
MAXIK = max(MAXI,MAXK)
! I1 and Xi1s must also be able to hold largest st block
! and (i,j,kstring) array

LSCR3 = max(MXIJST,MAXIK*MXTSOB,MXSTBL0)

call mma_allocate(I1,LSCR3,Label='I1')
call mma_allocate(I2,LSCR3,Label='I2')
call mma_allocate(I3,MAXIK*MXTSOB,Label='I3')
call mma_allocate(I4,MAXIK*MXTSOB,Label='I4')
call mma_allocate(XI1S,LSCR3,Label='XI1S')
call mma_allocate(XI2S,LSCR3,Label='XI2S')
call mma_allocate(XI3S,MAXIK*MXTSOB,Label='XI3S')
call mma_allocate(XI4S,MAXIK*MXTSOB,Label='XI4S')
LSCR2 = max(MXCJ,MXCIJA,MXCIJB,MXCIJAB)

! Merethe Interface
MOCAA = 0
if (MOCAA /= 0) then
  ! These blocks will also be used for single excitations so,
  ! To be removed   Largest block of SX excitations
  MAXE = max(NAEL,NBEL)
  MAXEl3 = min(MAXE,MXTSOB)
  MXSXST = (MXTSOB+1)*MAXEL3
  MXSXBL = MXSXST*MXSTBL0
# ifdef _DEBUGPRINT_
  write(u6,*) ' MXSXST,MXSXBL = ',MXSXST,MXSXBL
# endif
  LSCR2 = max(4*MXSXBL,LSCR2)
end if

LSCR12 = max(LSCR1,2*LSCR2)

if (IDIAG == 2) then
  ! PICO diagonalizer uses block KVEC3, use this as scratch block
  write(u6,*) 'Unchartered territory!'
  call Abend()
  !pC2 => VEC3 ! this is not clear yet.
  if (2*LSCR2 > LSCR1) then
    call mma_allocate(CJRES,LSCR2,Label='CJRES')
    pCJRES => CJRES
    call mma_allocate(SIRES,LSCR2,Label='SIRES')
    pSIRES => SIRES
  !else
  !  pCJRES => VEC3            ! ditto
  !  pSIRES => VEC3(1+LSCR2:)  ! ditto
  end if
else
  call mma_allocate(C2,LSCR12,Label='C2')
  pC2 => C2
  pCJRES => C2
  pSIRES => C2(1+LSCR2:)
end if

! Symmetry handling symmetry allowed/forbidden

! Out KSBLTP [nIrrep]

call ZBLTP(ISSM,nIrrep,IDC,SBLTP,SVST)
call ZBLTP(ICSM,nIrrep,IDC,CBLTP,SVST)
! 10 OOS arrays
nOOS = NOCTPA*NOCTPB*nIrrep
call mma_allocate(OOS,nOOS,10,Label='OOS')

iiCOPY = 1
! Transform C vector from CSF to SD basis
if (NOCSF == 0) then
  call CSDTVC_MCLR_1(C,HC,DTOC,CNSM(kic(1))%ICTS,icsm,iiCOPY)
  HC(1:NSDET) = Zero
end if

! Transform from combination scaling to determinant scaling
if ((IDC /= 1) .and. (ICISTR == 1)) &
  call SCDTC2_MCLR(C,ICSM,CBLTP,nIrrep,NOCTPA,NOCTPB,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,CIOIO,IDC,2,IDUMMY)

if (I12 == 2) then
  IDOH2 = 1
else
  IDOH2 = 0
end if

if (ICISTR == 1) then
  LLUC = 0
  LLUHC = 0
else
  LLUC = LUC
  LLUHC = LUHC
end if

if (ICISTR == 1) then
  call RASSG4(C,HC,CB,SB,pC2,CIOIO,SIOIO,ICSM,ISSM,CBLTP,SBLTP,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,NAEL,IATP,NBEL,IBTP,NOCTPA,NOCTPB, &
              nIrrep,NTSOB,IBTSOB,ITSOB,MAXK,MAXI,LSCR1,LSCR1,INSCR,pCJRES,pSIRES,Str(IATP)%EL1,Str(IATP)%EL3,Str(IBTP)%EL1, &
              Str(IBTP)%EL3,IDC,OOS(:,1),OOS(:,2),OOS(:,3),OOS(:,4),OOS(:,5),OOS(:,6),OOS(:,7),OOS(:,8),OOS(:,9),OOS(:,10),I1, &
              XI1S,I2,XI2S,I3,XI3S,I4,XI4S,IDOH2,SVST,PSSIGN,LLUC,LLUHC,IST,pCJRES,pSIRES,NOPARt,TimeDep)

else
  call SysHalt('sigmavec')

  ! IF WE USE DISK REPLACE THE FIRST VARIABLES LIKE THIS
  ! call RASSG4(C,HC,C,HC!,

end if

! Transform from combination scaling to determinant scaling
if ((IDC /= 1) .and. (ICISTR == 1)) &
  call SCDTC2_MCLR(HC,ISSM,SBLTP,nIrrep,NOCTPA,NOCTPB,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,CIOIO,IDC,1,IDUMMY)

! Transform HC vector from SD to CSF basis
if (NOCSF == 0) call CSDTVC_MCLR_2(C,HC,DTOC,CNSM(kic(2))%ICTS,ISSM,1)

! Eliminate local memory

call mma_deallocate(INSCR)
call mma_deallocate(SIOIO)
call mma_deallocate(CIOIO)
call mma_deallocate(SBLTP)
call mma_deallocate(CBLTP)

call mma_deallocate(I4)
call mma_deallocate(I3)
call mma_deallocate(I2)
call mma_deallocate(I1)
call mma_deallocate(XI4S)
call mma_deallocate(XI3S)
call mma_deallocate(XI2S)
call mma_deallocate(XI1S)
call mma_deallocate(OOS)
if (ICISTR == 1) then
  call mma_deallocate(CB)
  call mma_deallocate(SB)
end if
call mma_deallocate(SVST)
nullify(pC2,pCJRES,pCJRES,pSIRES)
call mma_deallocate(C2,safe='*')
call mma_deallocate(CJRES,safe='*')
call mma_deallocate(SIRES,safe='*')

end subroutine SigmaVec
