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

use Symmetry_Info, only: Mul
use Cholesky, only: nSym
use ChoMP2, only: DoDens, iFirstS, iT1am, LiT1am, LnOcc, LnT1am, nT1am, nVir
#ifdef _BUG_
use ChoMP2, only: iPQ_prod, LiPQprod, LnBatOrb, LnPQprod, nDel, nFro, nOcc, nPQ_prod
#endif
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: Vec(*)
real(kind=wp), intent(_OUT_) :: Srt(*)
integer(kind=iwp), intent(in) :: nVec, iSym, iBatch
integer(kind=iwp) :: iSyma, iSymi, iVec, kOff0, kOff1, kOff2, kOff3, lCp

if (.not. DoDens) then
  do iVec=1,nVec

    kOff0 = nT1am(iSym)*(iVec-1)
    kOff1 = LnT1am(iSym,iBatch)*(iVec-1)

    do iSymi=1,nSym

      iSyma = Mul(iSymi,iSym)
      if ((LnOcc(iSymi,iBatch) > 0) .and. (nVir(iSyma) > 0)) then
        lCp = nVir(iSyma)*LnOcc(iSymi,iBatch)
        kOff2 = kOff0+iT1am(iSyma,iSymi)+nVir(iSyma)*(iFirstS(iSymi,iBatch)-1)
        kOff3 = kOff1+LiT1am(iSyma,iSymi,iBatch)
        Srt(kOff3+1:kOff3+lCp) = Vec(kOff2+1:kOff2+lCp)
      end if

    end do

  end do
else
# ifdef _BUG_
  ! Special sorting for Mp2-density calculations where all integrals
  ! are used (for pure energy calculations only integrals of type
  ! (occ,vir|occ,vir) are needed).
  do iVec=1,nVec

    kOff0 = nPQ_prod(iSym)*(iVec-1)
    kOff1 = LnPQprod(iSym,iBatch)*(iVec-1)

    do iSymI=1,nSym

      iSymA = Mul(iSymI,iSym)
      if ((LnBatOrb(iSymI,iBatch) > 0) .and. (nFro(iSymA)+nOcc(iSymA)+nVir(iSymA)+nDel(iSymA) > 0)) then
        lCp = (nFro(iSymA)+nOcc(iSymA)+nVir(iSymA)+nDel(iSymA))*LnBatOrb(iSymI,iBatch)
        kOff2 = kOff0+iPQ_prod(iSymA,iSymI)+(nFro(iSymA)+nOcc(iSymA)+nVir(iSymA)+nDel(iSymA))*(iFirstS(iSymi,iBatch)-1)
        kOff3 = kOff1+LiPQprod(iSyma,iSymi,iBatch)
        Srt(kOff3+1:kOff3+lCp) = Vec(kOff2+1:kOff2+lCp)
      end if
    end do
  end do
# else
  ! BUG: iPQ_prod was never initialized
  call WarningMessage(2,'Sorry, but there is a bug in ChoMP2_Srt')
  call Abend()
# endif
end if

end subroutine ChoMP2_Srt
