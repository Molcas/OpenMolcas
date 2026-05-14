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

function ChoMP2_Setup_MemChk(LnT1am,LnPQprod,NumVec,nFrac,nSym,nBatch,Mem)
!
! Thomas Bondo Pedersen, Nov. 2004.
!
! Purpose: Check memory availability.

use ChoMP2, only: Laplace, SOS_mp2
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: ChoMP2_Setup_MemChk
integer(kind=iwp), intent(in) :: nSym, nBatch, LnT1am(nSym,nBatch), LnPQprod(nSym,nBatch), NumVec(nSym), nFrac(nSym), Mem
integer(kind=iwp) :: iBatch, iSym, jBatch, LiPQRSprod(8), LiT2am(8), LnPQRSprod, LnT2am, Nai, NumV
real(kind=wp) :: xDiff, xDim, xInt, xLeft, xMem, xNeed

if (Mem < 1) then
  ChoMP2_Setup_MemChk = .false.
  return
else
  ChoMP2_Setup_MemChk = .true.
end if

if (Laplace .and. SOS_MP2) then
  xMem = real(mem,kind=wp)
  xNeed = Zero
  do iBatch=1,nBatch
    do iSym=1,nSym
      Nai = LnT1am(iSym,iBatch)
      if ((Nai > 0) .and. (NumVec(iSym) > 0)) xNeed = max(xNeed,real(Nai,kind=wp)*real(NumVec(iSym),kind=wp))
    end do
  end do
  xLeft = xMem-xNeed
  if (xLeft < Zero) then
    ChoMP2_Setup_MemChk = .false.
    return
  end if
else
  do iSym=1,nSym
    if (nFrac(iSym) < 1) then
      ChoMP2_Setup_MemChk = .false.
      return
    end if
  end do
  xMem = real(mem,kind=wp)
  do jBatch=1,nBatch
    do iBatch=1,jBatch
      call ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)
      if (.false.) call ChoMP2_Energy_GetPQInd(LnPQRSprod,LiPQRSprod,iBatch,jBatch)
      if (.false.) then
        xInt = real(LnPQRSprod,kind=wp)
      else
        xInt = real(LnT2am,kind=wp)
      end if
      xLeft = xMem-xInt
      if ((xInt < One) .or. (xLeft < One)) then
        ChoMP2_Setup_MemChk = .false.
        return
      end if
      do iSym=1,nSym
        if (iBatch == jBatch) then
          if (.false.) then
            xDim = real(LnPQprod(iSym,iBatch),kind=wp)
          else
            xDim = real(LnT1am(iSym,iBatch),kind=wp)
          end if
        else
          if (.false.) then
            xDim = real(LnPQprod(iSym,iBatch))+real(LnPQprod(iSym,iBatch),kind=wp)
          else
            xDim = real(LnT1am(iSym,iBatch))+real(LnT1am(iSym,jBatch),kind=wp)
          end if
        end if
        if (nFrac(iSym) > NumVec(iSym)) then
          NumV = min(NumVec(iSym),1)
        else
          NumV = NumVec(iSym)/nFrac(iSym)
        end if
        xNeed = xDim*real(NumV,kind=wp)
        xDiff = xLeft-xNeed
        if (xDiff < One) then
          ChoMP2_Setup_MemChk = .false.
          return
        end if
      end do
    end do
  end do
end if

end function ChoMP2_Setup_MemChk
