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

subroutine Cho_GnVc_GetInt(xInt,lInt,nVecRS,iVecRS,ListSp,mSym,mPass,mmShl,iPass1,NumPass,NumSP)
!
! Purpose: compute integrals for NumPass integral passes starting at
!          pass iPass1.

use Cholesky, only: InfVec, IndRSh, nnShl, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lInt, mSym, mPass, nVecRS(mSym,mPass), iVecRS(mSym,mPass), mmShl, iPass1, NumPass
real(kind=wp), intent(inout) :: xInt(lInt)
integer(kind=iwp), intent(out) :: ListSP(mmShl), NumSP
integer(kind=iwp) :: iAB, iPass, iPass2, iShAB, iSP, iSym, iV, iV1, iV2, jShAB, lSewInt
integer(kind=iwp), allocatable :: SPTmp(:)
character(len=*), parameter :: SecNam = 'Cho_GnVc_GetInt'
integer(kind=iwp), external :: Cho_F2SP

! Initialization and input check.
! -------------------------------

if (NumPass < 1) then
  NumSP = 0
  return
end if

if (mSym /= nSym) call Cho_Quit('Input error [1] in '//SecNam,103)

if (iPass1 < 1) call Cho_Quit('Input error [2] in '//SecNam,103)

iPass2 = iPass1+NumPass-1
if (iPass2 > mPass) call Cho_Quit('Input error [3] in '//SecNam,103)

if (mmShl < nnShl) call Cho_Quit('Input error [4] in '//SecNam,103)

! Set up list of shell pairs to compute.
! --------------------------------------

call mma_allocate(SPTmp,nnShl,Label='SPTmp')
SPTmp(:) = 0
NumSP = 0
do iPass=iPass1,iPass2
  do iSym=1,nSym
    iV1 = iVecRS(iSym,iPass)
    iV2 = iV1+nVecRS(iSym,iPass)-1
    do iV=iV1,iV2
      iAB = InfVec(iV,1,iSym) ! addr in 1st reduced set
      jShAB = IndRsh(iAB) ! shell pair (full)
      iShAB = Cho_F2SP(jShAB) ! reduced shell pair
      if (iShAB > 0) then
        if (SPTmp(iShAB) == 0) then ! register SP
          SPTmp(iShAB) = 1
          NumSP = NumSP+1
          ListSP(NumSP) = iShAB
        end if
      else
        call Cho_Quit('SP not found in reduced list!',103)
      end if
    end do
  end do
end do
call mma_deallocate(SPTmp)

! Set memory used by Seward.
! --------------------------

call mma_maxDBLE(lSewInt)
call xSetMem_Ints(lSewInt)

! Loop through shell pair list.
! -----------------------------

do iSP=1,NumSP

  ! Get shell pair index.
  ! ---------------------

  iShAB = ListSP(iSP)

  ! Compute integral distribution (**|A B).
  ! ---------------------------------------

  call Cho_MCA_CalcInt_3(xInt,lInt,iShAB)

end do

! Deallocation.
! -------------

call xRlsMem_Ints()

end subroutine Cho_GnVc_GetInt
