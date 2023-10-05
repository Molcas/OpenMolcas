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
! Copyright (C) 2008, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2_Energy_GetPQInd(LnPQRSprod,LiPQRSprod,iBatch,jBatch)
!
! Jonas Bostrom, june 2008
!
! Purpose: setup (pq|rs) index arrays for calculating mp2_densities.

use Index_Functions, only: nTri_Elem
use Cholesky, only: nSym
use ChoMP2, only: ChoAlg, LnPQprod
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: LnPQRSprod, LiPQRSprod(8)
integer(kind=iwp), intent(in) :: iBatch, jBatch
integer(kind=iwp) :: iSym
character(len=14) :: String
character(len=*), parameter :: SecNam = 'ChoMP2_Energy_GetPQInd'

if (iBatch == jBatch) then
  LnPQRSprod = 0
  if (ChoAlg == 1) then
    do iSym=1,nSym
      LiPQRSprod(iSym) = LnPQRSprod
      LnPQRSprod = LnPQRSprod+nTri_Elem(LnPQprod(iSym,iBatch))
    end do
  else
    write(String,'(A8,I6)') 'ChoAlg =',ChoAlg
    call SysAbendMsg(SecNam,'ChoAlg out-of-bounds error!',String)
  end if
else
  LnPQRSprod = 0
  do iSym=1,nSym
    LiPQRSprod(iSym) = LnPQRSprod
    LnPQRSprod = LnPQRSprod+LnPQprod(iSym,iBatch)*LnPQprod(iSym,jBatch)
  end do
end if

end subroutine ChoMP2_Energy_GetPQInd
