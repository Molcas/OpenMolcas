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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine Freeze(TrMat,nTrMat,OneHam,mBT)
!***********************************************************************
!                                                                      *
!     purpose: Modify transformation matrix such atomic orbitals we    *
!              want to freeze are put in the first position            *
!                                                                      *
!     input:                                                           *
!       TrMat   : transformation matrix of length nTrMat               *
!                                                                      *
!     output:                                                          *
!       TrMat   : modified transformation matrix                       *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!***********************************************************************

use InfSCF, only: nSym, nBas, nFro, nOrb, nBT
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
#include "Molcas.fh"
integer nTrMat, mBT
real*8 TrMat(nTrMat), OneHam(mBT)
! Define local variables
integer iStart, iSym, nOr, iFro, iStrt, iSta, iBas, i, j, k, iSwap, iCMO, Ind1, Ind2
real*8 OEMin
real*8, dimension(:), allocatable :: Temp
integer MapBas(MxBas,MxSym)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Allocate temporary spaces
call mma_allocate(Temp,nBT,Label='Temp')

call dcopy_(nBT,OneHam,1,Temp,1)

! Form an array saying which atomic orbitals are frozen
iStart = 0
do iSym=1,nSym
  nOr = nOrb(iSym)
  do iFro=1,nFro(iSym)
    OEMin = 1.0d+06
    iStrt = iStart
    iSta = 0
    do iBas=1,nOr
      iStrt = iStrt+iBas
      if (Temp(iStrt) < OEMin) then
        MapBas(iFro,iSym) = iBas
        OEMin = Temp(iStrt)
        iSta = iStrt
      end if
    end do
    if (iSta /= 0) Temp(iSta) = -Temp(iSta)
  end do
  iStart = iStart+nOr*(nOr+1)/2
  do i=1,nFro(iSym)-1
    k = i
    do j=i+1,nFro(iSym)
      if (MapBas(j,iSym) < MapBas(k,iSym)) k = j
    end do
    if (k /= i) then
      iSwap = MapBas(k,iSym)
      MapBas(k,iSym) = MapBas(i,iSym)
      MapBas(i,iSym) = iSwap
    end if
  end do
end do

! Move frozen atomic orbitals to the first positions
iCMO = 1
do iSym=1,nSym
  do iFro=1,nFro(iSym)
    ind1 = iCMO+(iFro-1)*nBas(iSym)
    ind2 = iCMO+(MapBas(iFro,iSym)-1)*nBas(iSym)
    call DSwap_(nBas(iSym),TrMat(ind1),1,TrMat(ind2),1)
  end do
  iCMO = iCMO+nBas(iSym)*nOrb(iSym)
end do

! Deallocate memory
call mma_deallocate(Temp)

end subroutine Freeze
