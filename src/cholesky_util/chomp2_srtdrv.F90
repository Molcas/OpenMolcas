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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_SrtDrv(irc,DelOrig)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: presort Cholesky vectors according to batch structure in
!          MP2 program.
!
! DelOrig: input : flag for deleting original vector files.
!          output: flag to tell that at least 1 symmetry block has
!                  in fact been deleted.

use Cholesky, only: nSym, NumCho
use ChoMP2, only: DecoMP2, LnT1am, lUnit, lUnit_F, nBatch, nMP2Vec, nT1am
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp), intent(inout) :: DelOrig
integer(kind=iwp) :: iAdr, iBat, iBatch, iClos, iOpt, iSym, iTyp, iVec1, kChoMO, kSort, lChoMO, LnT1amx, lSort, lTot, lWrk, &
                     MinMem, nBat, nSrtVec, NumV, NumVec
real(kind=wp), allocatable :: Wrk(:)
character(len=*), parameter :: SecNam = 'ChoMP2_SrtDrv'

irc = 0
if (nBatch < 1) return

! Allocate available memory.
! --------------------------

call mma_maxDBLE(lWrk)
call mma_allocate(Wrk,lWrk,Label='Wrk')

! Set vector type (i.e., transformed vectors or vectors from (ai|bj)
! decomposition. Decide whether original files should be deleted.
! ------------------------------------------------------------------

if (DecoMP2) then
  iTyp = 2
else
  iTyp = 1
end if

if (DelOrig) then
  iClos = 3
else
  iClos = 2
end if
DelOrig = .false.

! Start symmetry loop.
! --------------------

do iSym=1,nSym

  ! Set number of vectors.
  ! ----------------------

  if (iTyp == 1) then
    nSrtVec = NumCho(iSym)
  else if (iTyp == 2) then
    nSrtVec = nMP2Vec(iSym)
  else
    irc = -1
    call Finish_this()
    return
  end if

  if ((nT1am(iSym) > 0) .and. (nSrtVec > 0)) then

    ! Set up vector batch.
    ! --------------------

    LnT1amx = 0
    do iBatch=1,nBatch
      LnT1amx = max(LnT1amx,LnT1am(iSym,iBatch))
    end do

    MinMem = nT1am(iSym)+LnT1amx
    NumVec = min(lWrk/MinMem,nSrtVec)
    if (NumVec < 1) then
      irc = 1
      call Finish_this()
      return
    else
      nBat = (nSrtVec-1)/NumVec+1
    end if

    ! Open full vector file.
    ! ----------------------

    call ChoMP2_OpenF(1,iTyp,iSym)

    ! Start batch loop.
    ! -----------------

    do iBat=1,nBat

      if (iBat == nBat) then
        NumV = nSrtVec-NumVec*(nBat-1)
      else
        NumV = NumVec
      end if
      iVec1 = NumVec*(iBat-1)+1

      ! Read batch of vectors.
      ! ----------------------

      lChoMO = nT1am(iSym)*NumV

      kChoMO = 1

      iOpt = 2
      lTot = lChoMO
      iAdr = nT1am(iSym)*(iVec1-1)+1
      call ddaFile(lUnit_F(iSym,iTyp),iOpt,Wrk(kChoMO),lChoMO,iAdr)

      ! Sort and write to disk.
      ! -----------------------

      kSort = 1+lChoMO
      lSort = lWrk-lChoMO
      do iBatch=1,nBatch
        lTot = LnT1am(iSym,iBatch)*NumV
        if (lSort < lTot) call SysAbendMsg(SecNam,'sort batch error','[0]')
        call ChoMP2_Srt(Wrk(kChoMO),Wrk(kSort),NumV,iSym,iBatch)
        call ChoMP2_OpenB(1,iSym,iBatch)
        iOpt = 1
        iAdr = LnT1am(iSym,iBatch)*(iVec1-1)+1
        call ddaFile(lUnit(iSym,iBatch),iOpt,Wrk(kSort),lTot,iAdr)
        call ChoMP2_OpenB(2,iSym,iBatch)
      end do

    end do

    ! Close (and possibly delete) full vector file.
    ! ---------------------------------------------

    call ChoMP2_OpenF(iClos,iTyp,iSym)
    DelOrig = iClos == 3

  end if

end do

call Finish_this()

contains

subroutine Finish_this()

  call mma_deallocate(Wrk)

end subroutine Finish_this

end subroutine ChoMP2_SrtDrv
