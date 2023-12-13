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

use Slapaf_Info, only: Analytic_Hessian, BMx, BSet, Degen, dqInt, dqInt_Aux, Gx, Gx0, HSet, iInt, iOptC, KtB, Lbl, lOld, MaxItr, &
                       NAC, nBVec, nFix, Numerical, qInt, Smmtrc
use Kriging_Mod, only: nSet
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nsAtom, nDim, nIter, mTR
real(kind=wp), intent(in) :: Coor(3,nsAtom)
integer(kind=iwp), intent(out) :: nQQ
integer(kind=iwp) :: i, iAtom, iInter, ix, ixyz, j, nRowH
logical(kind=iwp) :: Proc_dB
real(kind=wp), allocatable :: Degen2(:)

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

! iOptC(8) = constrained optimization
Proc_dB = HSet .and. (.not. lOld) .and. (Analytic_Hessian .or. Numerical .or. btest(iOptC,8))
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
  if (nSet > 1) call Force(nFix,Gx0(:,:,nIter),nsAtom,nQQ,BMx,nIter,dqInt_Aux(:,:,1),Lbl,Degen)
  if (nSet > 2) call Force(nFix,NAC(:,:,nIter),nsAtom,nQQ,BMx,nIter,dqInt_Aux(:,:,2),Lbl,Degen)
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
    !KtB(1:nDim,iInter) = KtB(1:nDim,iInter)/Sqrt(Degen2(1:nDim))
    KtB(1:nDim,iInter) = KtB(1:nDim,iInter)/Degen2(1:nDim)
  end do
  call mma_deallocate(Degen2)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine BMtrx_User_Defined
