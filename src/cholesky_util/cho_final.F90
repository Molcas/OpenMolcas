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

subroutine CHO_FINAL(WriteBookmarks)
!
! Purpose: Cholesky finalizations.

use Cholesky, only: BkmThr, BkmVec, Cho_AdrVec, Cho_Reord, CHOINICHECK, iSOShl, nBasT, nCol_BkmThr, nCol_BkmVec, nRow_BkmThr, &
                    nRow_BkmVec, nSym, ThrCom
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: WriteBookmarks
integer(kind=iwp) :: CHOISINI, IREO, l, NUMV(8)
integer(kind=iwp), allocatable :: BkmDim(:), iScratch(:)
real(kind=wp), allocatable :: Scratch(:)

! Write NUMCHO array, shell indices, and threshold to runfile.
! ------------------------------------------------------------

call CHO_P_GETGV(NUMV,NSYM)
call PUT_IARRAY('NUMCHO',NUMV,NSYM)
call PUT_IARRAY('iSOShl',ISOSHL,NBAST)
call PUT_DSCALAR('Cholesky Threshold',THRCOM)
#ifdef _DEBUGPRINT_
! This is needed in order for bookmark tests in cho_x_init to work
if (WriteBookmarks) call Put_lScalar('1C-CD',Cho_1Center)
#endif

! Write bookmarks to runfile.
! First, transpose array.
! ---------------------------

if (WriteBookmarks) then
  call mma_allocate(BkmDim,4,Label='BkmDim')
  BkmDim(1) = nCol_BkmVec
  BkmDim(2) = nRow_BkmVec
  BkmDim(3) = nCol_BkmThr
  BkmDim(4) = nRow_BkmThr
  call Put_iArray('Cholesky BkmDim',BkmDim,size(BkmDim))
  call mma_deallocate(BkmDim)
  if ((nRow_BkmVec > 0) .and. (nCol_BkmVec > 0) .and. (nRow_BkmThr > 0) .and. (nCol_BkmThr > 0)) then
    l = nRow_BkmVec*nCol_BkmVec
    call mma_allocate(iScratch,l,Label='iScratch')
    call iTrnsps(nSym,nCol_BkmVec,BkmVec,iScratch)
    call Put_iArray('Cholesky BkmVec',iScratch,l)
    call mma_deallocate(iScratch)
    call mma_deallocate(BkmVec)
    nRow_BkmVec = 0
    nCol_BkmVec = 0
    l = nRow_BkmThr*nCol_BkmThr
    call mma_allocate(Scratch,l,Label='Scratch')
    call Trnsps(nSym,nCol_BkmThr,BkmThr,Scratch)
    call Put_dArray('Cholesky BkmThr',Scratch,l)
    call mma_deallocate(Scratch)
    call mma_deallocate(BkmThr)
    nRow_BkmThr = 0
    nCol_BkmThr = 0
  end if
end if
if (allocated(BkmVec)) then
  call mma_deallocate(BkmVec)
  nRow_BkmVec = 0
  nCol_BkmVec = 0
end if
if (allocated(BkmThr)) then
  call mma_deallocate(BkmThr)
  nRow_BkmThr = 0
  nCol_BkmThr = 0
end if

! Write vector file address mode to runfile.
! ------------------------------------------

call PUT_ISCALAR('ChoVec Address',CHO_ADRVEC)

! Write reorder mark to runfile.
! ------------------------------

if (CHO_REORD) then
  IREO = 1
else
  IREO = 0
end if
call PUT_ISCALAR('Cholesky Reorder',IREO)

! Set initialization integer flag to "not set".
! ---------------------------------------------

CHOISINI = CHOINICHECK+1
call PUT_ISCALAR('ChoIni',CHOISINI)

end subroutine CHO_FINAL
