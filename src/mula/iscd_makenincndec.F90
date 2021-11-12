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

subroutine ISCD_MakenIncDec(n_max,nOrd,nOsc,lNMAT,lNINC,lNDEC,lBatch,nBatch,leftBatch,nIndex,Graph2,nMat,nInc,nDec)

use stdalloc, only: mma_allocate, mma_deallocate

implicit real*8(a-h,o-z)
integer nMat(nOsc,lBatch), nInc(nOsc,lBatch), nDec(nOsc,lBatch)
integer Graph2(n_max+1,n_max+1,nOsc)
integer n_max, nOrd, lBatch, nBatch, leftBatch
integer lNMAT, lNINC, lNDEC
#include "io_mula.fh"
integer nIndex(3,0:maxMax_n)
integer, allocatable :: iVecD(:), iVecI(:)

!GGt -------------------------------------------------------------------
!write(u6,*)
!write(u6,*) 'CGGt[ISCD_Mk_nIncDec] Infos:                   '
!write(u6,*) '     nMat(',nOsc,',',lBatch,')'
!write(u6,*) '     n_max,nOrd,nOsc==',n_max,nOrd,nOsc
!write(u6,*) '     lBatch,nBatch,leftBatch==',lBatch,nBatch,leftBatch
!write(u6,*) '----------------------------------------------'
!write(u6,*) '  The nIndex file:    '
!do i=1,nBatch+1
!  write(u6,*) i,': ',nIndex(1,i)
!end do
!  write(u6,*) '----------------------------------------------'
!call XFlush(u6)
!GGt -------------------------------------------------------------------

! Initialize

call mma_allocate(iVecI,nOsc,label='iVecI')
call mma_allocate(iVecD,nOsc,label='iVecD')

! Macrocycle iBatch
! Reading nMat

iIndex = 0
jIndex = 0
do iBatch=1,nBatch
  do ii=1,lBatch
    nInc(:,ii) = -1
    nDec(:,ii) = -1
  end do
  kIndex = nIndex(1,iBatch)
  !write(u6,*)'          iBatch=',iBatch,'  kIndex=',kIndex
  call iDaFile(lNMAT,2,nMat,nOsc*lBatch,kIndex)
  !GGt -----------------------------------------------------------------
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,': ',(nMat(k,i+1),k=1,nOsc)
  !end do
  !write(u6,*) '----------------------------------------------'
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
  do ii=1,lBatch

    ! Create nInc.

    iVecI(:) = nMat(:,ii)
    do j=1,nOsc
      iVecI(j) = iVecI(j)+1
      nInc(j,ii) = iDetnr(iVecI,Graph2,nOsc,n_max)
      iVecI(j) = iVecI(j)-1
    end do

    ! Create nDec.

    do j=1,nOsc
      if (nMat(j,ii) /= 0) then
        iVecD(:) = nMat(:,ii)
        iVecD(j) = iVecD(j)-1
        nDec(j,ii) = iDetnr(iVecD,Graph2,nosc,n_max)
        do iv=1,nOsc
          iVecD(iv) = iVecD(j)+1
        end do
      else
        nDec(j,ii) = -1
      end if
    end do

  end do

  nIndex(2,iBatch) = iIndex
  call iDaFile(lNINC,1,nInc,nOsc*lBatch,iIndex)
  nIndex(3,iBatch) = jIndex
  call iDaFile(lNDEC,1,nDec,nOsc*lBatch,jIndex)
  !GGt -----------------------------------------------------------------
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,':I',(nInc(k,i+1),k=1,nOsc)
  !  write(u6,*) i+(iBatch-1)*lBatch,':D',(nDec(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
  !  write(u6,*) '            nInc Written at ',nIndex(2,iBatch)
  !  write(u6,*) '            nDec Written at ',nIndex(3,iBatch)
  !  write(u6,*) '----------------------------------------------'
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------

end do

if (leftBatch > 0) then
  kIndex = nIndex(1,nBatch+1)
  !write(u6,*) '          nBatch+1',nBatch+1,'  kIndex=',kIndex
  call iDaFile(lNMAT,2,nMat,nOsc*lBatch,kIndex)
  !GGt -----------------------------------------------------------------
  !write(u6,*) '  --------- last Batch'
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,': ',(nMat(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------

  !write(u6,*) 'CGGt nBatch*lBatch,nOrd==',nBatch*lBatch,nOrd
  !call XFlush(u6)
  do ii=1,lBatch
    nInc(:,ii) = -1
    nDec(:,ii) = -1
  end do
  do iOrd=nBatch*lBatch,nOrd
    ii = 1+iOrd-nBatch*lBatch

    ! Create nInc.

    iVecI(:) = nMat(:,ii)
    do j=1,nOsc
      iVecI(j) = iVecI(j)+1
      nInc(j,ii) = iDetnr(iVecI,Graph2,nOsc,n_max)
      iVecI(j) = iVecI(j)-1
    end do

    ! Create nDec.

    do j=1,nOsc
      if (nMat(j,ii) /= 0) then
        iVecD(:) = nMat(:,ii)
        iVecD(j) = iVecD(j)-1
        nDec(j,ii) = iDetnr(iVecD,Graph2,nosc,n_max)
        do iv=1,nOsc
          iVecD(iv) = iVecD(j)+1
        end do
      else
        nDec(j,ii) = -1
      end if
    end do

  end do

  nIndex(2,nBatch+1) = iIndex
  call iDaFile(lNINC,1,nInc,nOsc*lBatch,iIndex)
  nIndex(3,nBatch+1) = jIndex
  call iDaFile(lNDEC,1,nDec,nOsc*lBatch,jIndex)
  !GGt -----------------------------------------------------------------
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,':I',(nInc(k,i+1),k=1,nOsc)
  !  write(u6,*) i+(iBatch-1)*lBatch,':D',(nDec(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
  !write(u6,*) '            nInc Written at ',nIndex(2,nBatch+1)
  !write(u6,*) '            nDec Written at ',nIndex(3,nBatch+1)
  !write(u6,*) '----------------------------------------------'
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
end if

!write(u6,*) '----------------------------------------------'
!call XFlush(u6)

call mma_deallocate(iVecI)
call mma_deallocate(iVecD)

return

end subroutine ISCD_MakenIncDec
!####
subroutine ISCD_ReloadNMAT(lnTabDim,nOrd,nOsc,lNMAT0,lNMAT,lBatch,nBatch,leftBatch,nIndex,nTabDim,nMat0,nMat)

implicit real*8(a-h,o-z)
integer nMat(nOsc,lBatch)
integer nMat0(nOsc), nTabDim(0:lnTabDim)
integer lnTabDim, nOrd, nOsc, lNMAT0, lNMAT
integer lBatch, nBatch, leftBatch
#include "io_mula.fh"
integer nIndex(3,0:maxMax_n)

! Initialize

!GGt -------------------------------------------------------------------
!write(u6,*)
!write(u6,*) 'CGGt[ISCD_ReloadNMAT] Infos:                   '
!write(u6,*) '     nMat(',nOsc,',',lBatch,')'
!write(u6,*) '     lnTabDim,nOrd,nOsc==',lnTabDim,nOrd,nOsc
!write(u6,*) '     lBatch,nBatch,leftBatch==',lBatch,nBatch,leftBatch
!write(u6,*) '     lnTabDim+1=',lnTabDim+1,':'
!do i=0,lnTabDim
!  write(u6,*) i,' read at ',nTabDim(i)
!  iIndex0 = nTabDim(i)
!  call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
!  write(u6,*) i,' read at',nTabDim(i),'  M:',(nMat0(j),j=1,nOsc)
!end do
!write(u6,*) '----------------------------------------------'
!call XFlush(u6)
!GGt -------------------------------------------------------------------
jIndex = 0
rewind(lNMAT0)
do iBatch=1,nBatch
  do ii=0,lBatch-1
    iOrd = ii+(iBatch-1)*lBatch
    iIndex0 = nTabDim(iOrd)
    call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
    do iOsc=1,nOsc
      nMat(iOsc,ii+1) = nMat0(iOsc)
    end do
  end do
  nIndex(1,iBatch) = jIndex
  call iDaFile(lNMAT,1,nMat,nOsc*lBatch,jIndex)
  !GGt -----------------------------------------------------------------
  !write(u6,*)'  --------- iBatch =',iBatch
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,':M',(nMat(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
end do
if (leftBatch > 0) then
  do iOrd=nBatch*lBatch,nOrd
    ii = iOrd-nBatch*lBatch
    iIndex0 = nTabDim(iOrd)
    call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
    do iOsc=1,nOsc
      nMat(iOsc,ii+1) = nMat0(iOsc)
    end do
  end do
  nIndex(1,nBatch+1) = jIndex
  call iDaFile(lNMAT,1,nMat,nOsc*lBatch,jIndex)
  !GGt -----------------------------------------------------------------
  !write(u6,*) '  --------- last Batch'
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,':M',(nMat(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
end if
!GGt -------------------------------------------------------------------
!write(u6,*) '----------------------------------------------'
!write(u6,*) '  The nIndex file:    '
!do i=1,nBatch+1
! write(u6,*) i,': ',nIndex(1,i)
!end do
!write(u6,*) '----------------------------------------------'
!call XFlush(u6)
!GGt -------------------------------------------------------------------

return

end subroutine ISCD_ReloadNMAT
