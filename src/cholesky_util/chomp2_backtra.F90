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
! Copyright (C) 2007, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_BackTra(iTyp,COcc,CVir,BaseName_AO,DoDiag,Diag)
!
! Thomas Bondo Pedersen, Dec. 2007.
!
! Purpose: Backtransform MO vectors to AO basis.
!          The result vectors are stored in lower triangular
!          storage.
!
! Note: do not call this routine directly; use ChoMP2_VectorMO2AO()
!       instead !!

use Symmetry_Info, only: Mul
use Cholesky, only: nBas, nSym
use ChoMP2, only: iAOVir, iT1am, iT1AOT, lUnit_F, nMP2Vec, nOcc, nT1am, nT1AOT, nVir
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iTyp
real(kind=wp), intent(in) :: COcc(*), CVir(*)
character(len=3), intent(in) :: BaseName_AO
logical(kind=iwp), intent(in) :: DoDiag
real(kind=wp), intent(_OUT_) :: Diag(*)
integer(kind=iwp) :: AlBe, iAB(8,8), iAdr, iOpt, iSym, iSyma, iSymAl, iSymb, iSymBe, iSymi, iVec, kAOVec, kCOcc, kCVir, kD, kDiag, &
                     kMOVec, kTemp, l_Buf, lU_AO, lVec, MaxInCore, na, nAB(8), nAB_Tot, nAl, ni, nVecInCore, nVecOnDisk
character(len=4) :: FullName_AO
real(kind=wp), allocatable :: AOVec(:), Buf(:), MOVec(:), Temp(:)
character(len=*), parameter :: SecNam = 'ChoMP2_BackTra'

! Set up index arrays.
! --------------------

iAB(:,:) = 0
nAB_Tot = 0
do iSym=1,nSym
  nAB(iSym) = 0
  do iSymb=1,nSym
    iSyma = Mul(iSymb,iSym)
    iAB(iSyma,iSymb) = nAB(iSym)
    nAB(iSym) = nAB(iSym)+nBas(iSyma)*nBas(iSymb)
  end do
  nAB_Tot = nAB_Tot+nAB(iSym)
end do

! Backtransform.
! --------------

if (DoDiag) Diag(1:nAB_Tot) = Zero

kDiag = 0
do iSym=1,nSym

  if ((nAB(iSym) >= 1) .and. (nMP2Vec(iSym) >= 1)) then

    iOpt = 1
    call ChoMP2_OpenF(iOpt,iTyp,iSym)
    write(FullName_AO,'(A3,I1)') BaseName_AO,iSym
    lU_AO = 7
    call daName_MF_WA(lU_AO,FullName_AO)

    call mma_allocate(AOVec,nAB(iSym),Label='AOVec')
    call mma_allocate(Temp,nT1AOT(iSym),Label='Temp')
    call mma_allocate(MOVec,nT1Am(iSym),Label='MOVec')

    call mma_maxDBLE(l_Buf)
    if (l_Buf < nAB(iSym)) then
      call SysAbendMsg(SecNam,'Insufficient memory!',' ')
    else
      call mma_allocate(Buf,l_Buf,Label='Buf')
    end if
    MaxInCore = min(l_Buf/nAB(iSym),nMP2Vec(iSym))

    nVecOnDisk = 0
    nVecInCore = 0
    do iVec=1,nMP2Vec(iSym)

      iOpt = 2
      iAdr = nT1Am(iSym)*(iVec-1)+1
      lVec = nT1Am(iSym)
      call ddaFile(lUnit_F(iSym,iTyp),iOpt,MOVec,lVec,iAdr)

      do iSymi=1,nSym
        iSyma = Mul(iSymi,iSym)
        iSymAl = iSyma
        kCVir = iAOVir(iSymAl,iSyma)+1
        kMOVec = 1+iT1Am(iSyma,iSymi)
        kTemp = 1+iT1AOT(iSymi,iSymAl)
        na = max(nVir(iSyma),1)
        nAl = max(nBas(iSymAl),1)
        ni = max(nOcc(iSymi),1)
        call DGEMM_('T','T',nOcc(iSymi),nBas(iSymAl),nVir(iSyma),One,MOVec(kMOVec),na,CVir(kCVir),nAl,Zero,Temp(kTemp),ni)
      end do

      do iSymBe=1,nSym
        iSymAl = Mul(iSymBe,iSym)
        iSymi = iSymBe
        kTemp = 1+iT1AOT(iSymi,iSymAl)
        kCOcc = iT1AOT(iSymi,iSymBe)+1
        kAOVec = 1+iAB(iSymAl,iSymBe)
        ni = max(nOcc(iSymi),1)
        nAl = max(nBas(iSymAl),1)
        call DGEMM_('T','N',nBas(iSymAl),nBas(iSymBe),nOcc(iSymi),One,Temp(kTemp),ni,COcc(kCOcc),ni,Zero,AOVec(kAOVec),nAl)
      end do

      if (DoDiag) then
        do AlBe=1,nAB(iSym)
          kD = kDiag+AlBe
          Diag(kD) = Diag(kD)+AOVec(AlBe)**2
        end do
      end if

      call dCopy_(nAB(iSym),AOVec,1,Buf(1+nVecInCore),MaxInCore)
      nVecInCore = nVecInCore+1

      if ((nVecInCore == MaxInCore) .or. (iVec == nMP2Vec(iSym))) then
        do AlBe=1,nAB(iSym)
          iOpt = 1
          iAdr = nMP2Vec(iSym)*(AlBe-1)+nVecOnDisk+1
          lVec = nVecInCore
          call ddaFile(lU_AO,iOpt,Buf(1+MaxInCore*(AlBe-1)),lVec,iAdr)
        end do
        nVecOnDisk = nVecOnDisk+nVecInCore
        nVecInCore = 0
      end if

    end do
#   ifdef _DEBUGPRINT_
    if ((nVecOnDisk /= nMP2Vec(iSym)) .or. (nVecInCore /= 0)) call SysAbendMsg(SecNam,'Logical bug detected!',' [1]')
#   endif

    call mma_deallocate(Buf)
    call mma_deallocate(MOVec)
    call mma_deallocate(Temp)
    call mma_deallocate(AOVec)

    call daClos(lU_AO)
    iOpt = 2
    call ChoMP2_OpenF(iOpt,iTyp,iSym)

  end if

  ! cycle loop point

  if (DoDiag) kDiag = kDiag+nAB(iSym)

end do

end subroutine ChoMP2_BackTra
