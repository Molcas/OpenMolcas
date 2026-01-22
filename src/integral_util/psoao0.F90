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
! Copyright (C) 1990,2015, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PSOAO0(nSO,MemPrm,MemMax,ipMem1,ipMem2,Mem1,Mem2,DoFock,nSD,iSD4)
!***********************************************************************
!                                                                      *
!  Object: to partition the SO and AO block. It will go to some length *
!          before it will start and break up the SO block. This will   *
!          reduce the total flop count. However, as we are breaking up *
!          the AO block this will affect the vectorization. Hence, at  *
!          some point it will actually be better to recompute the      *
!          primitives.                                                 *
!          Current strategy:                                           *
!          1. Start reducing the length of the primitives in the order *
!             lPrim,jPrim.                                             *
!          2. Reduce the size of the SO block by reducing the number of*
!             basis functions in the order lBas, jBas.                 *
!          3. Terminate run telling job max and min of additional      *
!             memory needed to perform the calculation.                *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified for scratch to FckAcc                           *
!                                                                      *
!             Modified for unified Work2 and Work3 block. Febr. 2015   *
!***********************************************************************

use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use lw_Info, only: lwInt, lwSqn, lwSyb
use Gateway_global, only: force_part_c, force_part_p
use RICD_Info, only: Cholesky, Do_RI
use k2_arrays, only: DoGrad_
use Symmetry_Info, only: nIrrep
use Breit, only: nComp
use Molcas, only: lCache
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSO, MemPrm, MemMax, ipMem1, nSD
integer(kind=iwp), intent(out) :: ipMem2, Mem1, Mem2
logical(kind=iwp), intent(in) :: DoFock
integer(kind=iwp), intent(inout) :: iSD4(0:nSD,4)
integer(kind=iwp) :: iBas, iBsInc, iCmp, iFact, IncVec, iPrim, iPrInc, jBas, jBsInc, jCmp, jPrim, jPrInc, kBas, kBsInc, kCmp, &
                     kPrim, kPrInc, kSOInt, la, lb, lBas, lBsInc, lc, lCmp, ld, lPack, lPrim, lPrInc, lSize, mab, mabcd, mabMax, &
                     mabMin, mcdMax, mcdMin, Mem0, MemAux, MemCon, MemFck, MemPck, MemPr, MemSp1, mijkl, na1a, na1b, na2a, na2b, &
                     na3a, na3b, nab, nabcd, nCache_, ncd, ne, nf, nijkl, nVec1, nVec2
logical(kind=iwp) :: Fail, QiBas, QjBas, QjPrim, QkBas, QlBas, QlPrim

la = iSD4(1,1)
lb = iSD4(1,2)
lc = iSD4(1,3)
ld = iSD4(1,4)

iCmp = iSD4(2,1)
jCmp = iSD4(2,2)
kCmp = iSD4(2,3)
lCmp = iSD4(2,4)

iBas = iSD4(3,1)
jBas = iSD4(3,2)
kBas = iSD4(3,3)
lBas = iSD4(3,4)

iPrim = iSD4(5,1)
jPrim = iSD4(5,2)
kPrim = iSD4(5,3)
lPrim = iSD4(5,4)

nab = iCmp*jCmp
ncd = kCmp*lCmp
mabMin = nTri3_Elem1(max(la,lb)-1)
mabMax = nTri3_Elem1(la+lb)-1
ne = (mabMax-mabMin+1)
mcdMin = nTri3_Elem1(max(lc,ld)-1)
mcdMax = nTri3_Elem1(lc+ld)-1
nf = (mcdMax-mcdMin+1)
! e0|f0 size
mabcd = ne*nf*nComp
! ab|cd size
nabcd = nab*ncd*nComp

! For test purposed we can force partitioning for both the
! contracted and primitive block

if (force_part_c) then
  iBsInc = (iBas+1)/2
  jBsInc = (jBas+1)/2
  kBsInc = (kBas+1)/2
  lBsInc = (lBas+1)/2
else
  iBsInc = iBas
  jBsInc = jBas
  kBsInc = kBas
  lBsInc = lBas
end if
if (force_part_p) then
  jPrInc = (jPrim+1)/2
  lPrInc = (lPrim+1)/2
else
  jPrInc = jPrim
  lPrInc = lPrim
end if
iPrInc = iPrim
kPrInc = kPrim

iFact = 1
if (.not. (Cholesky .or. Do_RI)) iFact = 4+3

do

  ! We can partition the blocks on all contracted indices
  ! and on the second and fourth primitive indices.
  ! We put priority in keeping an as large fraction of the contracted
  ! block as much preserved to limit the recompution of primitive
  ! integrals.

  QjPrim = .false.
  QlPrim = .true.
  QiBas = .false.
  QjBas = .false.
  QkBas = .false.
  QlBas = .false.
  Mem0 = MemMax

  mijkl = iPrInc*jPrInc*kPrInc*lPrInc
  nijkl = iBsInc*jBsInc*kBsInc*lBsInc

  ! Work1
  ! Memory for SO block. If petite list format is used there
  ! will be no SO block.

  kSOInt = nSO*nijkl
  Mem1 = iFact*kSOInt
  if (Mem1 == 0) Mem1 = 1
  if (nIrrep == 1) Mem1 = 1+(iFact-1)*nabcd*nijkl
  if (Mem1+1 > Mem0) then
    QjPrim = .false.
    QlPrim = .false.
    QiBas = .false.
    QjBas = .false.
    QkBas = .false.
    QlBas = .true.
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                Fail)
    if (Fail) then
      call WarningMessage(2,' Allocation failed for Work1')
      write(u6,*) Mem0,Mem1
      write(u6,*) iPrInc,iBsInc,kPrInc,kBsInc,jPrInc,jBsInc,lPrInc,lBsInc
      call Abend()
    end if
    cycle
  end if
  Mem0 = Mem0-Mem1-1

  ! Work2
  ! MemPr  : Scratch for Rys

  ! Memory for the Rys-Gauss procedure. This includes memory for
  ! all the intermediate arrays AND the {e0|f0} block.

  MemPr = MemPrm*mijkl

  ! MemAux : Auxiliary memory for partial contraction.

  ! If the primitive block is not full we need to accumulate primitive
  ! contributions to the contracted block. This will require a
  ! permanant pice of memory, during the computation of the primitive
  ! sublocks, where the incomplete block of contracted integrals are
  ! stored (Work4) - mabcd*nijkl.

  if ((jPrInc /= jPrim) .or. (lPrInc /= lPrim)) then
    MemAux = max(mabcd,nabcd)*nijkl
  else
    MemAux = 0
  end if
  !write(u6,*) 'MemAux:',MemAux

  ! MemCon : Scratch for the contraction step

  ! We need memory to the largest of the {e0|f0} or (e0|f0) block.
  ! Note that this correspond to that the A3 block in the second
  ! halftransformation is set to zero. We also assume that the
  ! integrals might be transformed to real spherical harmonic before
  ! the contraction step. This only happens for mijkl>=nijl
  !
  ! The routine which does the contraction does this in two
  ! half contraction steps with an intermediate which contains
  ! a subset of the quater or three quarter transformed integrals.
  ! Since we do not know the order of contraction at this time we will
  ! have to compute for the worst case scenario.
  !
  ! Contraction of the two first indices: iPrim->iBas & jPrim->jBas
  ! while the third and fourth indices are uncontracted.
  nCache_ = (3*lCache)/4-iPrim*iBas-jPrim*jBas
  ! Note that we do not know here the order of the contraction, hence
  ! we assume worst case scenario. The worst case is when lSize is
  ! the smallest.
  lSize = iPrInc*jPrInc+min(jPrInc*iBsInc,iPrInc*jBsInc)
  ! Length of the compound index.
  nVec1 = kPrInc*lPrInc*max(nabcd,mabcd)
  ! Compute the subsection length of the compound index.
  IncVec = min(max(1,nCache_/lSize),nVec1)
  ! Size of the initial integrals
  nA1a = max(nabcd,mabcd)*mijkl
  ! Size of the intemediate array
  nA2a = max(iPrim,jPrim)*IncVec
  ! Size of the half transformed integrals
  nA3a = iBsInc*jBsInc*nVec1
  !write(u6,*)
  !write(u6,*) 'IncVec,iPrim,jPrim:',IncVec,iPrim,jPrim
  !write(u6,*) 'nVec1,lSize=',nVec1,lSize
  !write(u6,*) 'nA1,nA2,nA3:',nA1a,nA2a,nA3a

  ! Contraction of the two last indices: kPrim->kBas & lPrim->lBas
  ! while the first and second indices are contracted.
  nCache_ = (3*lCache)/4-kPrim*kBas-lPrim*lBas
  lSize = kPrInc*lPrInc+min(lPrInc*kBsInc,kPrInc*lBsInc)
  nVec2 = iBsInc*jBsInc*max(nabcd,mabcd)
  IncVec = min(max(1,nCache_/lSize),nVec2)
  nA1b = nVec2*kPrInc*lPrInc
  nA2b = max(kPrim,lPrim)*IncVec
  nA3b = kBsInc*lBsInc*nVec2
  if (MemAux /= 0) nA3b = 0
  MemCon = max(nA1a,nA3b)+max(nA2a,nA2b)+max(nA3a,nA1b)
  !write(u6,*) 'IncVec,kPrim,lPrim:',IncVec,kPrim,lPrim
  !write(u6,*) 'nVec2,lSize=',nVec1,lSize
  !write(u6,*) 'nA1,nA2,nA3:',nA1b,nA2b,nA3b
  !write(u6,*) 'MemCon     :',MemCon

  ! Contraction of the two last indices: kPrim->kBas & lPrim->lBas
  ! while the first and second indices are uncontracted.
  nCache_ = (3*lCache)/4-kPrim*kBas-lPrim*lBas
  lSize = kPrInc*lPrInc+min(lPrInc*kBsInc,kPrInc*lBsInc)
  nVec1 = iPrInc*jPrInc*max(nabcd,mabcd)
  IncVec = min(max(1,nCache_/lSize),nVec1)
  nA1a = max(nabcd,mabcd)*mijkl
  nA2a = IncVec*max(kPrim,lPrim)
  nA3a = kBsInc*lBsInc*nVec1
  !write(u6,*) 'IncVec,kPrim,lPrim:',IncVec,kPrim,lPrim
  !write(u6,*) 'nVec1,lSize=',nVec1,lSize
  !write(u6,*) 'nA1,nA2,nA3:',nA1a,nA2a,nA3a

  ! Contraction of the two first indices: iPrim->iBas & jPrim->jBas
  ! while the third and fourth indices are contracted.
  nCache_ = (3*lCache)/4-iPrim*iBas-jPrim*jBas
  lSize = iPrInc*jPrInc+min(jPrInc*iBsInc,iPrInc*jBsInc)
  nVec2 = kBsInc*lBsInc*max(nabcd,mabcd)
  IncVec = min(max(1,nCache_/lSize),nVec2)
  nA1b = nVec2*iPrInc*jPrInc
  nA2b = IncVec*max(iPrim,jPrim)
  nA3b = iBsInc*jBsInc*nVec2
  if (MemAux /= 0) nA3b = 0
  MemCon = max(MemCon,max(nA1a,nA3b)+max(nA2a,nA2b)+max(nA3a,nA1b))
  !write(u6,*) 'IncVec,iPrim,jPrim:',IncVec,iPrim,jPrim
  !write(u6,*) 'nVec2,lSize=',nVec1,lSize
  !write(u6,*) 'nA1,nA2,nA3:',nA1b,nA2b,nA3b
  !write(u6,*) 'MemCon     :',MemCon

  ! MemSp1 : Scratch for the transformation from cartesian to
  !          spherical harmonics. This is executed inside the
  !          primitive block only if (1) the set of primitive
  !          integrals are fewer than the contracted functions
  !          and (2) there is no partitioning on the primitive
  !          index.

  if (Dograd_) then
    mab = nTri_Elem1(la)*nTri_Elem1(lb)
  else
    mab = nab
  end if
  MemSp1 = (max(mabcd,nabcd)+mab*nf*nComp)*nijkl

  ! MemFck : Scratch for FckAcc.

  ! Memory to manupulate the 1-particle densities in FckAcc.

  if (DoFock) then
    MemFck = 6*max(iBsInc*lBsInc,iBsInc*kBsInc,jBsInc*lBsInc,jBsInc*kBsInc,iBsInc*jBsInc,kBsInc*lBsInc)+nijkl*nabcd
  else
    MemFck = 0
  end if

  ! MemPck : Scratch for integral packing

  ! Memory for integral packing.

  if (.not. (Cholesky .or. Do_RI)) then
    MemPck = 2*nabcd*nijkl
  else
    MemPck = 0
  end if
  Mem2 = max(max(MemPr,MemCon,MemSp1)+MemAux,MemFck,MemPck)
  if (Mem2+1 > Mem0) then
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                Fail)
    if (Fail) then
      call WarningMessage(2,' Allocation failed for Work2')
      write(u6,*) Mem0
      write(u6,*) iPrInc,iBsInc,kPrInc,kBsInc,jPrInc,jBsInc,lPrInc,lBsInc
      call Abend()
    end if
  else
    exit
  end if
end do
Mem0 = Mem0-Mem2-1

ipMem2 = ipMem1+Mem1

! Auxiliary memory for integral packing

if (.not. (Cholesky .or. Do_RI)) then
  if (nIrrep == 1) then
    lPack = nabcd*nijkl
    lwInt = ipMem1
  else
    lPack = kSOInt
    lwInt = ipMem1+lPack
  end if
  lwSyB = lwInt+2*lPack
  lwSqN = lwSyB+2*lPack
else
  lPack = 0
  lwInt = 0
  lwSyB = 0
  lwSqN = 0
end if

iSD4(4,1) = iBsInc
iSD4(4,2) = jBsInc
iSD4(4,3) = kBsInc
iSD4(4,4) = lBsInc

iSD4(6,1) = iPrInc
iSD4(6,2) = jPrInc
iSD4(6,3) = kPrInc
iSD4(6,4) = lPrInc

end subroutine PSOAO0
