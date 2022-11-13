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
! Copyright (C) 2013, Roland Lindh                                     *
!               2015, Ignacio Fdez. Galvan (split from gencxctl)       *
!***********************************************************************

subroutine get_drdq(drdq,mInt,nLambda,mLambda,Iter)
!***********************************************************************
! subroutine to get the dr/dq vectors for the constraints as given in  *
! the 'UDC' file.                                                      *
!***********************************************************************

use Slapaf_Info, only: BMx, Degen
use Slapaf_Parameters, only: iRow_c, Curvilinear
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mInt, nLambda, Iter
real(kind=wp), intent(inout) :: drdq(mInt,nLambda)
integer(kind=iwp), intent(out) :: mLambda
integer(kind=iwp) :: n3, nBV, i, iLambda, iOff, iOff2, iAtom, ixyz
real(kind=wp) :: RR
logical(kind=iwp) :: lWrite
integer(kind=iwp), allocatable :: iFlip(:)
real(kind=wp), allocatable :: BVc(:), dBVc(:), BMx_t(:,:), Val(:), Val0(:), cInt(:), cInt0(:), Mult(:), dBMx(:)
character(len=8), allocatable :: Lbl(:)
real(kind=wp), external :: DDot_

lWrite = .false.
n3 = size(Degen)

if (nLambda /= 0) then
  nBV = iRow_c-nLambda-1
  call mma_allocate(BMx_t,n3,nLambda,Label='BMx_t')

  call mma_allocate(BVc,n3*nBV,Label='BVc')
  call mma_allocate(dBVc,nBV*n3**2,Label='dBVc')
  call mma_allocate(Val,nBV,Label='Val')
  call mma_allocate(Val0,nBV,Label='Val0')
  Val0(:) = Zero
  call mma_allocate(cInt,nLambda,Label='cInt')
  call mma_allocate(cInt0,nLambda,Label='cInt0')
  call mma_allocate(Mult,nBV**2,Label='Mult')
  call mma_allocate(dBMx,nLambda*n3**2,Label='dBMx')
  call mma_allocate(iFlip,nBV,Label='iFlip')
  call mma_allocate(Lbl,mInt,Label='Lbl')

  call DefInt2(BVc,dBVc,nBV,BMx_t,nLambda,size(Degen,2),iRow_c,Val,cInt,cInt0,Lbl,lWrite,Mult,dBMx,Val0,Iter,iFlip)

  call mma_deallocate(Lbl)
  call mma_deallocate(iFlip)
  call mma_deallocate(dBMx)
  call mma_deallocate(Mult)
  call mma_deallocate(cInt0)
  call mma_deallocate(cInt)
  call mma_deallocate(Val0)
  call mma_deallocate(Val)
  call mma_deallocate(dBVc)
  call mma_deallocate(BVc)

# ifdef _DEBUGPRINT_
  call RecPrt('BMx_t',' ',BMx_t,n3,nLambda)
# endif

  ! Assemble dr/dq: Solve  B dr/dq = dr/dx

  drdq(:,:) = Zero

  ! Temporary fix of the dC/dx vector which always
  ! is propted up with the full degeneracy factor.

  if (.not. Curvilinear) then
    do iLambda=1,nLambda
      do i=1,n3
        iAtom = (i+2)/3
        ixyz = i-(iAtom-1)*3
        BMx_t(i,iLambda) = BMx_t(i,iLambda)/Degen(ixyz,iAtom)
      end do
    end do
  end if

  call Eq_Solver('N',n3,mInt,nLambda,BMx,Curvilinear,Degen,BMx_t,drdq)
# ifdef _DEBUGPRINT_
  call RecPrt('drdq',' ',drdq,mInt,nLambda)
# endif

  call mma_deallocate(BMx_t)
end if

! Double check that we don't have any null vectors

iOff = 1
iOff2 = 1
mLambda = nLambda
do iLambda=1,nLambda
  RR = sqrt(DDot_(mInt,drdq(1,iOff),1,drdq(1,iOff),1))
  if (RR < 1.0e-12_wp) then
    write(u6,*) 'Warning: constraint ',iLambda,' has a null vector, I''ll remove it!'
    mLambda = mLambda-1
  else
    if (iOff /= iOff2) call dCopy_(mInt,drdq(1,iOff),1,drdq(1,iOff2),1)
    iOff2 = iOff2+1
  end if
  iOff = iOff+1
end do
if (mLambda < nLambda) call FZero(drdq(1,mLambda+1),mInt*(nLambda-mLambda))

return

end subroutine get_drdq
