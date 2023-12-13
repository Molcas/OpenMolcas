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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine PSOAO0_h(nSO,nMemab,nMemcd,MemPrm,MemMax,iAnga,iCmpa,iBas,iBsInc,jBas,jBsInc,kBas,kBsInc,lBas,lBsInc,iPrim,iPrInc, &
                    jPrim,jPrInc,kPrim,kPrInc,lPrim,lPrInc,ipMem1,ipMem2,ipMem3,ipMem4,Mem1,Mem2,Mem3,Mem4,Mend)
!***********************************************************************
!                                                                      *
!  Object: to partion the SO and AO block. It will go to some length   *
!          before it will start and break up the SO block. This will   *
!          reduce the total flop count. However, as we are breaking up *
!          the AO block this will affect the vectorization. Hence, at  *
!          some point it will actually be better to recompute the      *
!          primitives.                                                 *
!          Current stratergy:                                          *
!          1. Start reducing the length of the primitives in the order *
!             lPrim,jPrim.                                             *
!          2. Reduce the size of the SO block by reducing the number of*
!             basis functions in the order lBas, jBas.                 *
!          3. Terminate run telling job max and min of additional      *
!             memory needed to perform the calculation.                *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!***********************************************************************

use Index_Functions, only: nTri3_Elem1, nTri_Elem1
use Symmetry_Info, only: nIrrep
use Definitions, only: iwp, u6, RtoI

implicit none
integer(kind=iwp), intent(in) :: nSO, nMemab, nMemcd, MemPrm, MemMax, iAnga(4), iCmpa(4), iBas, jBas, kBas, lBas, iPrim, jPrim, &
                                 kPrim, lPrim, ipMem1
integer(kind=iwp), intent(out) :: iBsInc, jBsInc, kBsInc, lBsInc, iPrInc, jPrInc, kPrInc, lPrInc, ipMem2, ipMem3, ipMem4, Mem1, &
                                  Mem2, Mem3, Mem4, Mend
#include "Molcas.fh"
#include "warnings.h"
integer(kind=iwp) :: iCmp, iFact, IncVec, jCmp, kCmp, kSOInt, la, lb, lc, lCmp, ld, lSize, mabcd, mabMax, mabMin, mcdMax, mcdMin, &
                     Mem0, MemAux, MemCon, MemPr, MemSp1, MemSp2, MemTr1, MemTr2, MemTr3, nA2, nA3, nCache_, nVec1, nVec2
logical(kind=iwp) :: Fail, QiBas, QjBas, QjPrim, QkBas, QlBas, QlPrim

!iQ = 1
la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
iCmp = iCmpa(1)
jCmp = iCmpa(2)
kCmp = iCmpa(3)
lCmp = iCmpa(4)
mabMin = nTri3_Elem1(max(la,lb)-1)
mabMax = nTri3_Elem1(la+lb)-1
mcdMin = nTri3_Elem1(max(lc,ld)-1)
mcdMax = nTri3_Elem1(lc+ld)-1
mabcd = (mabMax-mabMin+1)*(mcdMax-mcdMin+1)

iBsInc = iBas
jBsInc = jBas
kBsInc = kBas
lBsInc = lBas
iPrInc = iPrim
jPrInc = jPrim
kPrInc = kPrim
lPrInc = lPrim
iFact = 1
iFact = 4/RtoI+3

do
  QjPrim = .false.
  QlPrim = .true.
  QiBas = .false.
  QjBas = .false.
  QkBas = .false.
  QlBas = .false.
  Mem0 = MemMax

  ! Work1
  ! Memory for SO block. If petite list format is used there
  ! will be no SO block.

  kSOInt = nSO*iBsInc*jBsInc*kBsInc*lBsInc
  Mem1 = iFact*kSOInt
  if (Mem1 == 0) Mem1 = 1
  if (nIrrep == 1) Mem1 = 1+(iFact-1)*iCmp*jCmp*kCmp*lCmp*iBsInc*jBsInc*kBsInc*lBsInc
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
      write(u6,*) ' Allocation failed for Work1'
      write(u6,*) Mem0,Mem1
      write(u6,*) iPrInc,iBsInc,kPrInc,kBsInc,jPrInc,jBsInc,lPrInc,lBsInc
      call Quit(_RC_MEMORY_ERROR_)
    end if
    cycle
  end if
  !write(u6,*) ' After Mem1',iPrInc,iBsInc,kPrInc,kBsInc,jPrInc,jBsInc,lPrInc,lBsInc
  Mem0 = Mem0-Mem1-1

  ! Work2
  ! MemPr  : Scratch for Rys
  ! MemCon : Scratch for the contraction step
  ! MemTr1 : Scratch for the 1st application of the transfer eqn.
  ! MemTr2 : Scratch for the 2nd application of the transfer eqn.

  ! Work4 (this is overlayed with Work2 and is placed at the top)
  ! MemAux : Auxiliary memory for partial contraction.

  ! Work5 (this is overlayed with Work2 and is placed at the top)

  MemPr = MemPrm*iPrInc*jPrInc*kPrInc*lPrInc
  MemCon = mabcd*max(iPrInc*jPrInc*kPrInc*lPrInc,iBsInc*jBsInc*kBsInc*lBsInc)
  if ((jPrInc /= jPrim) .or. (lPrInc /= lPrim)) then
    MemAux = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
  else
    MemAux = 0
  end if
  MemTr1 = (mabMax-mabMin+1)*nMemcd*iBsInc*jBsInc*kBsInc*lBsInc
  MemTr2 = kCmp*lCmp*nMemab*iBsInc*jBsInc*kBsInc*lBsInc
  Mem2 = max(MemPr+MemAux,MemCon+MemAux,MemTr1,MemTr2)
  if (Mem2+1 > Mem0) then
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                Fail)
    if (Fail) then
      write(u6,*) ' Allocation failed for Work2'
      write(u6,*) Mem0,Mem2,MemPr+MemAux,MemCon+MemAux,MemTr1,MemTr2
      write(u6,*) iPrInc,iBsInc,kPrInc,kBsInc,jPrInc,jBsInc,lPrInc,lBsInc
      call Quit(_RC_MEMORY_ERROR_)
    end if
    cycle
  end if
  if ((jPrInc /= jPrim) .or. (lPrInc /= lPrim)) then
    Mem4 = MemAux
  else
    Mem4 = Mem2
  end if
  !write(u6,*) ' After Mem2',jPrInc,jBsInc,lPrInc,lBsInc
  Mem0 = Mem0-Mem2-1

  ! Work3
  ! MemCon : Scratch for the contraction step, and transpose after the contraction step
  ! MemSp1 : Scratch for the 1st transformation from cartesian to spherical harmonics.
  ! MemSp2 : Scratch for the 2nd transformation from cartesian to spherical harmonics.
  ! MemTr3 : Scratch for transpose in tnsctl

  nCache_ = (3*lCache)/4-iPrim*iBas-jPrim*jBas
  lSize = iPrInc*jPrInc+max(jPrInc*iBsInc,iPrInc*jBsInc)
  nVec1 = kPrInc*lPrInc*mabcd
  IncVec = min(max(1,nCache_/lSize),nVec1)
  nA2 = max(jPrInc*iBsInc,iPrInc*jBsInc)*IncVec
  nA3 = iBsInc*jBsInc*nVec1
  MemCon = max(mabcd*iBsInc*jBsInc*kBsInc*lBsInc,nA2+nA3)

  nCache_ = (3*lCache)/4-kPrim*kBas-lPrim*lBas
  lSize = kPrInc*lPrInc+max(lPrInc*kBsInc,kPrInc*lBsInc)
  nVec2 = iBsInc*jBsInc*mabcd
  IncVec = min(max(1,nCache_/lSize),nVec2)
  nA2 = max(lPrInc*kBsInc,kPrInc*lBsInc)*IncVec
  !nA3 = kBsInc*lBsInc*nVec2
  MemCon = max(MemCon,nA3+nA2)

  !*********************************************************************

  nCache_ = (3*lCache)/4-kPrim*kBas-lPrim*lBas
  lSize = kPrInc*lPrInc+max(lPrInc*kBsInc,kPrInc*lBsInc)
  nVec1 = iPrInc*jPrInc*mabcd
  IncVec = min(max(1,nCache_/lSize),nVec1)
  nA2 = IncVec*max(lPrInc*kBsInc,kPrInc*lBsInc)
  nA3 = kBsInc*lBsInc*nVec1
  MemCon = max(MemCon,nA3+nA2)

  nCache_ = (3*lCache)/4-iPrim*iBas-jPrim*jBas
  lSize = iPrInc*jPrInc+max(jPrInc*iBsInc,iPrInc*jBsInc)
  nVec2 = kBsInc*lBsInc*mabcd
  IncVec = min(max(1,nCache_/lSize),nVec2)
  nA2 = IncVec*max(jPrInc*iBsInc,iPrInc*jBsInc)
  !nA3 = iBsInc*jBsInc*nVec2
  MemCon = max(MemCon,nA3+nA2)

  MemSp1 = (mabMax-mabMin+1)*lCmp*nTri_Elem1(lc)*iBsInc*jBsInc*kBsInc*lBsInc
  MemSp2 = lCmp*kCmp*jCmp*nTri_Elem1(la)*iBsInc*jBsInc*kBsInc*lBsInc
  MemTr3 = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
  Mem3 = max(MemCon,MemSp1,MemSp2,MemTr3)
  if (Mem3+1 <= Mem0) exit
  call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
              Fail)
  if (Fail) then
    write(u6,*) ' Allocation failed for Work3'
    write(u6,*) Mem0,Mem3,MemCon,MemSp1,MemSp2
    write(u6,*) iPrInc,iBsInc,kPrInc,kBsInc,jPrInc,jBsInc,lPrInc,lBsInc
    call Quit(_RC_MEMORY_ERROR_)
  end if
end do
Mem0 = Mem0-Mem3-1

ipMem2 = ipMem1+Mem1
ipMem3 = ipMem2+Mem2
ipMem4 = ipMem2+Mem2-Mem4
Mend = 0

return

end subroutine PSOAO0_h
