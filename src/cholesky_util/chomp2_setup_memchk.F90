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

logical function ChoMP2_Setup_MemChk(LnT1am,LnPQprod,NumVec,nFrac,nSym,nBatch,Mem)
!
! Thomas Bondo Pedersen, Nov. 2004.
!
! Purpose: Check memory availability.

implicit real*8(a-h,o-z)
#include "chomp2_cfg.fh"
integer LnT1am(nSym,nBatch)
integer LnPQprod(nSym,nBatch)
integer NumVec(nSym), nFrac(nSym)
integer LnT2am, LiT2am(8)
integer LnPQRSprod, LiPQRSprod(8)
logical Accepted

if (Mem < 1) then
  Accepted = .false.
  Go To 1 ! exit
else
  Accepted = .true.
end if

if (Laplace .and. SOS_MP2) then
  xMem = dble(mem)
  xNeed = 0.0d0
  do iBatch=1,nBatch
    do iSym=1,nSym
      Nai = LnT1am(iSym,iBatch)
      if ((Nai > 0) .and. (NumVec(iSym) > 0)) xNeed = max(xNeed,dble(Nai)*dble(NumVec(iSym)))
    end do
  end do
  xLeft = xMem-xNeed
  if (xLeft < 1.0d0) then
    Accepted = .false.
    Go To 1 ! exit
  end if
else
  do iSym=1,nSym
    if (nFrac(iSym) < 1) then
      Accepted = .false.
      Go To 1 ! exit
    end if
  end do
  xMem = dble(mem)
  do jBatch=1,nBatch
    do iBatch=1,jBatch
      call ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)
      if (.false.) call ChoMP2_Energy_GetPQInd(LnPQRSprod,LiPQRSprod,iBatch,jBatch)
      if (.false.) then
        xInt = dble(LnPQRSprod)
      else
        xInt = dble(LnT2am)
      end if
      xLeft = xMem-xInt
      if ((xInt < 1.0d0) .or. (xLeft < 1.0d0)) then
        Accepted = .false.
        Go To 1 ! exit
      end if
      do iSym=1,nSym
        if (iBatch == jBatch) then
          if (.false.) then
            xDim = dble(LnPQprod(iSym,iBatch))
          else
            xDim = dble(LnT1am(iSym,iBatch))
          end if
        else
          if (.false.) then
            xDim = dble(LnPQprod(iSym,iBatch))+dble(LnPQprod(iSym,iBatch))
          else
            xDim = dble(LnT1am(iSym,iBatch))+dble(LnT1am(iSym,jBatch))
          end if
        end if
        if (nFrac(iSym) > NumVec(iSym)) then
          NumV = min(NumVec(iSym),1)
        else
          NumV = NumVec(iSym)/nFrac(iSym)
        end if
        xNeed = xDim*dble(NumV)
        xDiff = xLeft-xNeed
        if (xDiff < 1.0d0) then
          Accepted = .false.
          Go To 1 ! exit
        end if
      end do
    end do
  end do
end if

1 ChoMP2_Setup_MemChk = Accepted

end function ChoMP2_Setup_MemChk
