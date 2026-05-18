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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine TRDNS1(IVEC,DPT1,NDPT1)
! Add to the transition density matrix DPT1,
!    DPT1(p,q) = Add <IVEC| E(p,q) |0>.
! where IVEC stands for the 1st-order perturbed CASPT2
! wave function stored as vector nr IVEC on LUSOLV.
! DPT1 is stored as symmetry-blocked array of square matrices.
! Each square matrix is actually lower block triangle, but is
! stored in full, including zero elements.

! Only cases A, C and D(Symm 1) contributes.

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
#endif
use fake_GA, only: GA_Arrays
use caspt2_module, only: nActel, nAsh, nAsh, nASup, nInDep, nIsh, nISup, nOrb, nSsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, NDPT1
real(kind=wp), intent(inout) :: DPT1(NDPT1)
integer(kind=iwp) :: IA, ICASE, ID, IDOFF, II, IMLTOP, ISYM, IT, ITTOT, IW, IWAI, IWAT, IWOFF, IWTI, lVec, NA, NAS, NI, NIS, NO, &
                     NS, NVEC, NWAI, NWAT, NWTI
real(kind=wp) :: FACT
real(kind=wp), allocatable :: WAI(:), WAT(:), WTI(:)
#ifdef _MOLCAS_MPP_
real(kind=wp), allocatable :: TMP(:)
#include "global.fh"
#include "mafdecls.fh"
#endif

! Transform to standard representation, covariant form.
call PTRTOC(1,IVEC,IVEC)

NWTI = 0
NWAI = 0
NWAT = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NA = NASH(ISYM)
  NS = NSSH(ISYM)
  NWTI = NWTI+NA*NI
  NWAI = NWAI+NS*NI
  NWAT = NWAT+NS*NA
end do

IMLTOP = 1
if (NWTI == 0) goto 110
call mma_allocate(WTI,NWTI,LABEL='WTI')
WTI(:) = Zero
ICASE = 1
IWOFF = 1
do ISYM=1,NSYM
  if (NINDEP(ISYM,ICASE) == 0) goto 100
  NIS = NISUP(ISYM,ICASE)
  NAS = NASUP(ISYM,ICASE)
  NVEC = NIS*NAS
  if (NVEC == 0) goto 100
  call RHS_ALLO(NAS,NIS,LVEC)
  call RHS_READ(NAS,NIS,LVEC,ICASE,ISYM,IVEC)
  FACT = One/real(max(1,NACTEL),kind=wp)
# ifdef _MOLCAS_MPP_
  if (IS_REAL_PAR()) then
    if (KING()) then
      call mma_allocate(TMP,NVEC,Label='TMP')
      call RHS_GET(NAS,NIS,LVEC,TMP)
      call SPEC1A(IMLTOP,FACT,ISYM,TMP,size(TMP),WTI(IWOFF),size(WTI(IWOFF:)))
      call mma_deallocate(TMP)
    end if
  else
# endif
    call SPEC1A(IMLTOP,FACT,ISYM,GA_Arrays(LVEC)%A,size(GA_Arrays(LVEC)%A),WTI(IWOFF),size(WTI(IWOFF:)))
# ifdef _MOLCAS_MPP_
  end if
# endif
  call RHS_FREE(LVEC)
  IWOFF = IWOFF+NASH(ISYM)*NISH(ISYM)
100 continue
end do
110 continue

if (NWAT == 0) goto 210
call mma_allocate(WAT,NWAT,Label='WAT')
WAT(:) = Zero
ICASE = 4
IWOFF = 1
do ISYM=1,NSYM
  if (NINDEP(ISYM,ICASE) == 0) goto 200
  NIS = NISUP(ISYM,ICASE)
  NAS = NASUP(ISYM,ICASE)
  NVEC = NIS*NAS
  if (NVEC == 0) goto 200
  if (NSSH(ISYM)*NASH(ISYM) == 0) goto 200
  call RHS_ALLO(NAS,NIS,LVEC)
  call RHS_READ(NAS,NIS,LVEC,ICASE,ISYM,IVEC)
  FACT = One/real(max(1,NACTEL),kind=wp)
# ifdef _MOLCAS_MPP_
  if (IS_REAL_PAR()) then
    if (KING()) then
      call mma_allocate(TMP,NVEC,LABEL='TMP')
      call RHS_GET(NAS,NIS,LVEC,TMP)
      call SPEC1C(IMLTOP,FACT,ISYM,TMP,size(TMP),WAT(IWOFF),size(WAT(IWOFF:)))
      call mma_deallocate(TMP)
    end if
  else
# endif
    call SPEC1C(IMLTOP,FACT,ISYM,GA_Arrays(LVEC)%A,size(GA_Arrays(LVEC)%A),WAT(IWOFF),size(WAT(IWOFF:)))
# ifdef _MOLCAS_MPP_
  end if
# endif
  call RHS_FREE(LVEC)
  IWOFF = IWOFF+NSSH(ISYM)*NASH(ISYM)
200 continue
end do
210 continue

if (NWAI == 0) goto 300
call mma_allocate(WAI,NWAI,Label='WAI')
WAI(:) = Zero
ICASE = 5
ISYM = 1
if (NINDEP(ISYM,ICASE) == 0) goto 300
NIS = NISUP(ISYM,ICASE)
NAS = NASUP(ISYM,ICASE)
NVEC = NIS*NAS
if (NVEC == 0) goto 300
call RHS_ALLO(NAS,NIS,LVEC)
call RHS_READ(NAS,NIS,LVEC,ICASE,ISYM,IVEC)
FACT = One/real(max(1,NACTEL),kind=wp)
#ifdef _MOLCAS_MPP_
if (IS_REAL_PAR()) then
  if (KING()) then
    call mma_allocate(TMP,NVEC,LABEL='TMP')
    call RHS_GET(NAS,NIS,LVEC,TMP)
    call SPEC1D(IMLTOP,FACT,TMP,NVEC,WAI,NWAI)
    call mma_deallocate(TMP)
  end if
else
#endif
  call SPEC1D(IMLTOP,FACT,GA_Arrays(LVEC)%A,size(GA_Arrays(LVEC)%A),WAI,nWAI)
#ifdef _MOLCAS_MPP_
end if
#endif
call RHS_FREE(LVEC)
300 continue

! Transform vectors back to eigenbasis of H0(diag).
call PTRTOSR(0,IVEC,IVEC)

if (NWTI > 0) call GADGOP(WTI,NWTI,'+')
if (NWAI > 0) call GADGOP(WAI,NWAI,'+')
if (NWAT > 0) call GADGOP(WAT,NWAT,'+')
! Put transition density elements in temporaries W into
! proper positions, as subdiagonal matrices in DPT1:
IDOFF = 0
IWTI = 0
IWAI = 0
IWAT = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NA = NASH(ISYM)
  NS = NSSH(ISYM)
  NO = NORB(ISYM)
  do IT=1,NA
    ITTOT = IT+NI
    do II=1,NI
      ID = IDOFF+ITTOT+NO*(II-1)
      IW = IWTI+IT+NA*(II-1)
      DPT1(ID) = DPT1(ID)+WTI(IW)
    end do
  end do
  do IT=1,NA
    ITTOT = IT+NI
    ID = IDOFF+NI+NA+NO*(ITTOT-1)
    IW = IWAT+IT-NA
    do IA=1,NS
      ID = ID+1
      IW = IW+NA
      DPT1(ID) = DPT1(ID)+WAT(IW)
    end do
  end do
  do II=1,NI
    ID = IDOFF+NI+NA+NO*(II-1)
    IW = IWAI+II-NI
    do IA=1,NS
      ID = ID+1
      IW = IW+NI
      DPT1(ID) = DPT1(ID)+WAI(IW)
    end do
  end do
  IWTI = IWTI+NA*NI
  IWAI = IWAI+NS*NI
  IWAT = IWAT+NS*NA
  IDOFF = IDOFF+NO**2
end do

if (NWTI > 0) call mma_deallocate(WTI)
if (NWAI > 0) call mma_deallocate(WAI)
if (NWAT > 0) call mma_deallocate(WAT)

end subroutine TRDNS1
