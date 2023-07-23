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

subroutine ChoMP2_Srt(Vec,Srt,nVec,iSym,iBatch)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: copy out subblock of vectors.

use ChoMP2, only: iFirstS, LiPQprod, LiT1am, LnBatOrb, LnOcc, LnPQprod, LnT1am
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Vec(*), Srt(*)
integer(kind=iwp) :: nVec, iSym, iBatch
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
integer(kind=iwp) :: iSyma, iSymi, iVec, kOff0, kOff1, kOff2, kOff3, lCp
! Statement function
integer(kind=iwp) :: MulD2h, i, j
MulD2h(i,j) = ieor(i-1,j-1)+1

if (.not. DoDens) then
  do iVec=1,nVec

    kOff0 = nT1am(iSym)*(iVec-1)+1
    kOff1 = LnT1am(iSym,iBatch)*(iVec-1)+1

    do iSymi=1,nSym

      iSyma = MulD2h(iSymi,iSym)
      if ((LnOcc(iSymi,iBatch) > 0) .and. (nVir(iSyma) > 0)) then
        lCp = nVir(iSyma)*LnOcc(iSymi,iBatch)
        kOff2 = kOff0+iT1am(iSyma,iSymi)+nVir(iSyma)*(iFirstS(iSymi,iBatch)-1)
        kOff3 = kOff1+LiT1am(iSyma,iSymi,iBatch)
        call dCopy_(lCp,Vec(kOff2),1,Srt(kOff3),1)
      end if

    end do

  end do
else
  ! Special sorting for Mp2-density calculations where all integrals
  ! are used (for pure energy calculations only integrals of type
  ! (occ,vir|occ,vir) are needed).
  do iVec=1,nVec

    kOff0 = nPQ_prod(iSym)*(iVec-1)+1
    kOff1 = LnPQprod(iSym,iBatch)*(iVec-1)+1

    do iSymI=1,nSym

      iSymA = MulD2h(iSymI,iSym)
      if ((LnBatOrb(iSymI,iBatch) > 0) .and. (nFro(iSymA)+nOcc(iSymA)+nVir(iSymA)+nDel(iSymA) > 0)) then
        lCp = (nFro(iSymA)+nOcc(iSymA)+nVir(iSymA)+nDel(iSymA))*LnBatOrb(iSymI,iBatch)
        kOff2 = kOff0+iPQ_prod(iSymA,iSymI)+(nFro(iSymA)+nOcc(iSymA)+nVir(iSymA)+nDel(iSymA))*(iFirstS(iSymi,iBatch)-1)
        kOff3 = kOff1+LiPQprod(iSyma,iSymi,iBatch)
        call dCopy_(lCp,Vec(kOff2),1,Srt(kOff3),1)
      end if
    end do
  end do
end if

end subroutine ChoMP2_Srt
