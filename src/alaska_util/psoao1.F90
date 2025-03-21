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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine PSOAO1(nSO,MemPrm,MemMax,iFnc,ipMem1,ipMem2,Mem1,Mem2,MemPSO,nSD,iSD4)
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
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to first order derivatives. January '92         *
!***********************************************************************

use PSO_Stuff, only: Gamma_On, lPSO, nGamma
use SOAO_Info, only: iAOtSO
use Gateway_global, only: force_part_c, force_part_p
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSO, MemPrm, MemMax, ipMem1, nSD
integer(kind=iwp), intent(out) :: iFnc(4), ipMem2, Mem1, Mem2, MemPSO
integer(kind=iwp), intent(inout) :: iSD4(0:nSD,4)
integer(kind=iwp) :: i1, iAO(4), iBas, iBsInc, iCmp, iCmpa(4), iFac, iiBas(4), IncVec, iPrim, iPrInc, iTmp1, j, jBas, jBsInc, &
                     jCmp, jPam, jPrim, jPrInc, kBas, kBsInc, kCmp, kPrim, kPrInc, kSOInt, la, lb, lBas, lBsInc, lc, lCmp, ld, &
                     lPrim, lPrInc, lSize, mabcd, Mem0, Mem3, MemAux, MemAux0, MemDeP, MemRys, MemScr, MemSph, MemTrn, nA2, nA3, &
                     nabcd, nCache, nFac, nPam(4,0:7), nTmp1, nTmp2, nVec1
logical(kind=iwp) :: Fail, QiBas, QjBas, QjPrim, QkBas, QlBas, QlPrim
integer(kind=iwp), external :: MemTra
#include "Molcas.fh"

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

iAO(:) = iSD4(7,:)
iCmpa(:) = iSD4(2,:)

mabcd = nTri_Elem1(la)*nTri_Elem1(lb)*nTri_Elem1(lc)*nTri_Elem1(ld)
nabcd = iCmp*jCmp*kCmp*lCmp

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

do
  QjPrim = .false.
  QlPrim = .true.
  QiBas = .false.
  QjBas = .false.
  QkBas = .false.
  QlBas = .false.
  Mem0 = MemMax

  ! *** Work1 ***

  ! Memory for 2nd order density matrix in SO basis.

  kSOInt = nSO*iBsInc*jBsInc*kBsInc*lBsInc
  Mem1 = kSOInt

  ! Allocate memory for MO to SO/AO transformation
  ! of the 2nd order density matrix for this shell quadruplet.

  if (lPSO) then
    iiBas(1) = iBsInc
    iiBas(2) = jBsInc
    iiBas(3) = kBsInc
    iiBas(4) = lBsInc
    nPam(:,:) = 0
    MemPSO = 1
    nTmp2 = 0

    do jPam=1,4
      iTmp1 = 0
      nTmp1 = 0
      do j=0,nIrrep-1
        do i1=1,iCmpa(jPam)
          if (iAOtSO(iAO(jPam)+i1,j) > 0) then
            nPam(jPam,j) = nPam(jPam,j)+iiBas(jPam)
            nTmp1 = nTmp1+iiBas(jPam)
            iTmp1 = iTmp1+1
          end if
        end do
      end do
      MemPSO = MemPSO*nTmp1
      nTmp2 = nTmp2+nTmp1
      iFnc(jPam) = iTmp1
    end do
    MemScr = MemTra(nPam)
    nFac = 4
    nTmp2 = nTmp2+4
  else
    MemScr = 0
    MemPSO = 0
    nFac = 0
    nTmp2 = 0
  end if
  MemAux0 = MemPSO+MemScr+nFac*S%nDim+nTmp2+4
  if (Mem1+1+MemAux0 > Mem0) then
    QjPrim = .false.
    QlPrim = .false.
    QiBas = .false.
    QjBas = .false.
    QkBas = .false.
    QlBas = .true.
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                Fail)
    if (Fail) then
      write(u6,*) ' Memory allocation failed for Work1'
      write(u6,'(2I3,L1,2I3,L1)') iBas,iBsInc,QiBas,kBas,kBsInc,QkBas
      write(u6,'(2I3,L1,2I3,L1)') jBas,jBsInc,QjBas,lBas,lBsInc,QlBas
      write(u6,'(2I3,L1,2I3,L1)') jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim
      write(u6,*) MemMax,Mem0,Mem1,MemAux0+1
      write(u6,*) MemPSO,MemScr,4*S%nDim,nTmp2+4
      call Abend()
    end if
    cycle
  end if
  ! Subtract one additional word (?)
  Mem0 = Mem0-Mem1-1

  ! *** Work2 and Work4 ***

  ! Memory for 2nd order density matrix in contracted basis (both
  ! cartesian and spherical harmonic) and in primitive basis.
  ! MemDeP: Target for desymmetrization
  ! MemTrn: Scratch and target for decontraction
  ! MemAux: Contracted 2nd order density matrix (if partial decon.)
  ! MemSph: transformation spherical harmonics to cartesian, source
  !         and target.
  ! nGamma: temporary storage of bins as read from file.

  if (Gamma_On) then
    iiBas(1) = iBas
    iiBas(2) = jBas
    iiBas(3) = kBas
    iiBas(4) = lBas
    nGamma = 1

    do jPam=1,4
      nTmp1 = 0
      do j=0,nIrrep-1
        do i1=1,iCmpa(jPam)
          if (iAOtSO(iAO(jPam)+i1,j) > 0) then
            nTmp1 = nTmp1+iiBas(jPam)
          end if
        end do
      end do
      nGamma = nGamma*nTmp1
    end do
  else
    nGamma = 0
  end if

  MemDeP = nabcd*iBsInc*jBsInc*kBsInc*lBsInc

  MemTrn = mabcd*max(iBsInc*jBsInc*kBsInc*lBsInc,iPrInc*jPrInc*kPrInc*lPrInc)
  ! If partial decontraction we need to keep the contracted 2nd
  ! order density matrix. (Work4)
  if ((jPrInc /= jPrim) .or. (lPrInc /= lPrim)) then
    MemAux = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
  else
    MemAux = 0
  end if
  MemSph = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
  Mem2 = max(MemTrn+MemAux,MemDeP,MemSph,nGamma+MemAux0)
  if (Mem2+1 > Mem0) then
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                Fail)
    if (Fail) then
      write(u6,*) ' Memory allocation failed for Work2'
      write(u6,'(2I3,L1,2I3,L1)') iBas,iBsInc,QiBas,kBas,kBsInc,QkBas
      write(u6,'(2I3,L1,2I3,L1)') jBas,jBsInc,QjBas,lBas,lBsInc,QlBas
      write(u6,'(2I3,L1,2I3,L1)') jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim
      write(u6,*) MemMax,Mem0,Mem1,Mem2
      call Abend()
    end if
    cycle
  end if
  ! Subtract one additional word (?)
  Mem0 = Mem0-Mem2-1

  ! *** Work3 and Work5 ***

  ! Scratch for decontraction and transformation to spherical gaussian.
  ! Working array for Rysg1.
  ! Scratch area for resolving degeneracies due to the double coset
  ! treatement of the symmetry.
  ! MemTrn: Scratch for decontraction
  ! MemRys: Scratch for calualation of primitive integral gradients.
  ! MemScr: Scratch for prescreening

  ! If mabcd /= 1 we need to do a transpose!
  iFac = 0
  if (mabcd /= 1) iFac = 1

  nCache = (3*lCache)/4-iBas*iPrim-jBas*jPrim
  lSize = iBsInc*jBsInc+max(iBsInc*jPrInc,iPrInc*jBsInc)
  nVec1 = kBsInc*lBsInc*mabcd
  IncVec = min(max(1,nCache/lSize),nVec1)
  nA2 = max(iBsInc*jPrInc,iPrInc*jBsInc)*IncVec
  nA3 = iPrInc*jPrInc*nVec1
  MemTrn = nA2+nA3

  nCache = (3*lCache)/4-kBas*kPrim-lBas*lPrim
  lSize = kBsInc*lBsInc+max(kBsInc*lPrInc,kPrInc*lBsInc)
  nVec1 = iPrInc*jPrInc*mabcd
  IncVec = min(max(1,nCache/lSize),nVec1)
  nA2 = max(kBsInc*lPrInc,kPrInc*lBsInc)*IncVec
  MemTrn = max(MemTrn,nA3+nA2,iFac*mabcd*iBsInc*jBsInc*kBsInc*lBsInc)
  !                                                                      *
  !***********************************************************************
  !                                                                      *
  nCache = (3*lCache)/4-kBas*kPrim-lBas*lPrim
  lSize = kBsInc*lBsInc+max(kBsInc*lPrInc,kPrInc*lBsInc)
  nVec1 = iBsInc*jBsInc*mabcd
  IncVec = min(max(1,nCache/lSize),nVec1)
  nA2 = max(kBsInc*lPrInc,kPrInc*lBsInc)*IncVec
  nA3 = kPrInc*lPrInc*nVec1
  MemTrn = max(MemTrn,nA3+nA2)

  MemTrn = max(MemTrn,nA2+nA3)
  nCache = (3*lCache)/4-iBas*iPrim-jBas*jPrim
  lSize = iBsInc*jBsInc+max(iBsInc*jPrInc,iPrInc*jBsInc)
  nVec1 = kPrInc*lPrInc*mabcd
  IncVec = min(max(1,nCache/lSize),nVec1)
  nA2 = max(iBsInc*jPrInc,iPrInc*jBsInc)*IncVec

  MemTrn = max(MemTrn,nA2+nA3,iFac*mabcd*iBsInc*jBsInc*kBsInc*lBsInc)

  MemRys = MemPrm*iPrInc*jPrInc*kPrInc*lPrInc
  MemScr = (2*mabcd+1)*iPrInc*jPrInc*kPrInc*lPrInc+iPrInc*jPrInc+kPrInc*lPrInc
  Mem3 = max(MemTrn,MemRys,MemScr)
  if (Mem3+1 > Mem0) then
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                Fail)
    if (Fail) then
      write(u6,*) ' Memory allocation failed for Work3'
      write(u6,'(2I3,L1,2I3,L1)') iBas,iBsInc,QiBas,kBas,kBsInc,QkBas
      write(u6,'(2I3,L1,2I3,L1)') jBas,jBsInc,QjBas,lBas,lBsInc,QlBas
      write(u6,'(2I3,L1,2I3,L1)') jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim
      write(u6,*) MemMax,Mem0,Mem1,Mem2,Mem3
      write(u6,*) MemTrn,MemRys,MemPrm
      call Abend()
    end if
  else
    exit
  end if
end do
! Subtract one additional word (?)
Mem0 = Mem0-Mem3-1

Mem2 = Mem2+Mem3

ipMem2 = ipMem1+Mem1

iSD4(4,1) = iBsInc
iSD4(4,2) = jBsInc
iSD4(4,3) = kBsInc
iSD4(4,4) = lBsInc

iSD4(6,1) = iPrInc
iSD4(6,2) = jPrInc
iSD4(6,3) = kPrInc
iSD4(6,4) = lPrInc

end subroutine PSOAO1
