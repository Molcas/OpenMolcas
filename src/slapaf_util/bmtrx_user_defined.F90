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

subroutine BMtrx_User_Defined(nsAtom,Coor,nDim,nIter,mTR,nQQ)

use Slapaf_Info, only: Gx, qInt, dqInt, KtB, BMx, Degen, Smmtrc, Lbl, Gx0, dqInt_Aux, NAC
use Slapaf_Parameters, only: iInt, nFix, nBVec, Analytic_Hessian, MaxItr, iOptC, BSet, HSet, lOld, Numerical
use Kriging_Mod, only: nSet

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
real*8 Coor(3,nsAtom)
logical Proc_dB
real*8, allocatable :: Degen2(:)
!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Section for user defined internal coordinates

call Rd_UDIC(iInt,nFix,nRowH) ! nRowH is not used!
nQQ = iInt+nFix

if (allocated(qInt)) then
  if (size(qInt,1) /= nQQ) then
    call mma_deallocate(qInt)
    call mma_deallocate(dqInt)
  end if
end if
if (.not. allocated(qInt)) then
  call mma_allocate(qInt,nQQ,MaxItr,Label='qInt')
  call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
  qInt(:,:) = Zero
  dqInt(:,:) = Zero
end if
if (allocated(dqInt_Aux)) then
  if (size(dqInt_Aux,1) /= nQQ) call mma_deallocate(dqInt_Aux)
end if
if ((.not. allocated(dqInt_Aux)) .and. (nSet > 1)) then
  call mma_allocate(dqInt_Aux,nQQ,MaxItr,nSet-1,Label='dqInt_Aux')
  dqInt_Aux(:,:,:) = Zero
end if
call mma_allocate(BMx,3*nsAtom,nQQ,Label='BMx')
BMx(:,:) = Zero

! Compute the B matrix in symmetry distinct basis and the
! internal coordinates.

! iOptC(256) = constrained optimization
Proc_dB = HSet .and. (.not. lOld) .and. (Analytic_Hessian .or. Numerical .or. (iand(iOptC,256) == 256))
! Compute and store dBQQ in the reference structure
if (Proc_dB) then
  ! Not implimented, sorry
end if

call DefInt(nBVec,BMx,nQQ,nsAtom,qInt(:,nIter),Lbl,Coor,nDim-mTR)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the gradient

if (BSet) then
  call Force(nFix,Gx(:,:,nIter),nsAtom,nQQ,BMx,nIter,dqInt,Lbl,Degen)
  if (nSet > 1) then
    call Force(nFix,Gx0(:,:,nIter),nsAtom,nQQ,BMx,nIter,dqInt_Aux(:,:,1),Lbl,Degen)
  end if
  if (nSet > 2) then
    call Force(nFix,NAC(:,:,nIter),nsAtom,nQQ,BMx,nIter,dqInt_Aux(:,:,2),Lbl,Degen)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (HSet .and. (.not. lOld) .and. BSet) then
  call mma_allocate(KtB,nDim,nQQ,Label='KtB')

  call mma_allocate(Degen2,nDim,Label='Degen2')
  i = 0
  do ix=1,3*nsAtom
    iAtom = (ix+2)/3
    ixyz = ix-(iAtom-1)*3
    if (Smmtrc(ixyz,iAtom)) then
      i = i+1
      Degen2(i) = Degen(ixyz,iAtom)
    end if
  end do

  do j=1,nQQ
    i = 0
    do ix=1,3*nsAtom
      iAtom = (ix+2)/3
      ixyz = ix-(iAtom-1)*3
      if (Smmtrc(ixyz,iAtom)) then
        i = i+1
        KtB(i,j) = BMx(ix,j)
      end if
    end do
  end do

  do iInter=1,nQQ
    do iDim=1,nDim
      !KtB(iDim,iInter) = KtB(iDim,iInter)/Sqrt(Degen2(iDim))
      KtB(iDim,iInter) = KtB(iDim,iInter)/Degen2(iDim)
    end do
  end do
  call mma_deallocate(Degen2)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine BMtrx_User_Defined
