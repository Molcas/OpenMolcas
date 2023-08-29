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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2g_Reord_R(Wrk,lWrk)
!
! Jonas Bostrom, Apr 2010
!
! Purpose: To reorder R-vectors so it is practical to access
!          one ia-piece at the time.

use Symmetry_Info, only: Mul
use Cholesky, only: nSym
use ChoMP2, only: AdrR1, AdrR2, iT1am, lUnit_F, LuRInv, nMoMo, nMP2Vec, nOcc, nVir
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lWrk
real(kind=wp), intent(out) :: Wrk(lWrk)
integer(kind=iwp) :: iA, iAdr, iAdr1, iAdr2, iBat, iClos, iI, ioffset1, iOffset2, iOpt, iSeed, iSym, iSymA, iSymI, iTypR, iVec, &
                     iVec1, iVecOV, kEndRia1, kRia1, kRia2, lRia, lTot, maxvalue, nBatR, NumVec, nVec
character(len=5) :: Fname
character(len=*), parameter :: SecNam = 'ChoMP2g_Reord_r'
integer(kind=iwp), external :: IsFreeUnit

iTypR = 2
iVecOV = 6
maxvalue = 1000

! Do not delete vectors
! ---------------------
iClos = 2

iSeed = 7
LuRInv(1) = IsFreeUnit(iSeed)
write(Fname,'(A4,I1)') 'TMPV',4
call DaName_MF_WA(LuRInv(1),Fname)

LuRInv(2) = IsFreeUnit(iSeed)
write(Fname,'(A4,I1)') 'TMPV',5
call DaName_MF_WA(LuRInv(2),Fname)

iAdr1 = 1
iAdr2 = 1
do iSymI=1,nSym
  do iSymA=1,nSym
    iSym = Mul(iSymA,iSymI)
    do iI=1,nOcc(iSymI)
      AdrR1(iSymA,iSymI,iI) = iAdr1
      iAdr1 = iAdr1+nVir(iSymA)*nMP2Vec(iSym)
    end do
    do iA=1,nVir(iSymA)
      AdrR2(iSymA,iSymI,iA) = iAdr2
      iAdr2 = iAdr2+nOcc(iSymI)*nMP2Vec(iSym)
    end do
  end do
end do

do iSym=1,nSym
  if (nMP2Vec(iSym) == 0) cycle
  nVec = min(maxvalue,nMP2Vec(iSym))
  if (nVec < 1) call SysAbendMsg(SecNam,'Insufficient memory','[1]')
  nBatR = (nMP2Vec(iSym)-1)/nVec+1

  ! Allocate memory for Ria-vectors
  ! -------------------------------

  lRia = nMoMo(iSym,iVecOV)*nVec
  kRia1 = 1
  kEndRia1 = kRia1+lRia

  kRia2 = kEndRia1

  ! Open Cholesky amplitude vectors
  ! -------------------------------
  call ChoMP2_OpenF(1,iTypR,iSym)

  do iBat=1,nBatR
    if (iBat == nBatR) then
      NumVec = nMP2Vec(iSym)-nVec*(nBatR-1)
    else
      NumVec = nVec
    end if
    iVec = nVec*(iBat-1)+1

    ! Read Amplitude vectors
    ! ----------------------
    iOpt = 2
    lTot = nMoMo(iSym,iVecOV)*NumVec
    iAdr = nMoMo(iSym,iVecOV)*(iVec-1)+1
    call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia1),lTot,iAdr)

    do iVec1=1,NumVec
      do iSymI=1,nSym
        iSymA = Mul(iSymI,iSym)
        do iI=1,nOcc(iSymI)
          ioffset1 = (iI-1)*nVir(iSymA)+iT1am(iSymA,iSymI)+(iVec1-1)*nMoMo(iSym,iVecOV)
          iOffset2 = (iVec1-1)*nVir(iSymA)+(iI-1)*NumVec*nVir(iSymA)+iT1am(iSymA,iSymI)*NumVec
          Wrk(kRia2+iOffset2:kRia2+iOffset2+nVir(iSymA)-1) = Wrk(kRia1+iOffset1:kRia1+iOffset1+nVir(iSymA)-1)
        end do
      end do
    end do

    ! Put the reordered vectors on disk
    iOpt = 1
    do iSymI=1,nSym
      iSymA = Mul(iSymI,iSym)
      do iI=1,nOcc(iSymI)
        lTot = nVir(iSymA)*NumVec
        iAdr = AdrR1(iSymA,iSymI,iI)+(iVec-1)*nVir(iSymA)
        iOffset2 = (iI-1)*NumVec*nVir(iSymA)+iT1am(iSymA,iSymI)*NumVec
        call dDaFile(LuRInv(1),iOpt,Wrk(kRia2+iOffSet2),lTot,iAdr)
      end do
    end do

  end do     !iBat
  do iBat=1,nBatR
    if (iBat == nBatR) then
      NumVec = nMP2Vec(iSym)-nVec*(nBatR-1)
    else
      NumVec = nVec
    end if
    iVec = nVec*(iBat-1)+1

    ! Read Amplitude vectors
    ! ----------------------
    iOpt = 2
    lTot = nMoMo(iSym,iVecOV)*NumVec
    iAdr = nMoMo(iSym,iVecOV)*(iVec-1)+1
    call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia1),lTot,iAdr)

    do iSymI=1,nSym
      iSymA = Mul(iSymI,iSym)
      do iI=1,nOcc(iSymI)
        do iA=1,nVir(iSymA)

          iOffSet1 = iA-1+(iI-1)*nVir(iSymA)+iT1am(iSymA,iSymI)
          iOffSet2 = iI-1+(iA-1)*NumVec*nOcc(iSymI)+iT1am(iSymA,iSymI)*NumVec
          call dCopy_(NumVec,Wrk(kRia1+iOffset1),nMoMo(iSym,iVecOV),Wrk(kRia2+iOffset2),nOcc(iSymI))
        end do
      end do
    end do

    ! Put the reordered vectors on disk
    iOpt = 1
    do iSymI=1,nSym
      iSymA = Mul(iSymI,iSym)
      do iA=1,nVir(iSymA)
        lTot = nOcc(iSymI)*NumVec
        iAdr = AdrR2(iSymA,iSymI,iA)+(iVec-1)*nOcc(iSymI)
        iOffset2 = (iA-1)*NumVec*nOcc(iSymI)+iT1am(iSymA,iSymI)*NumVec
        call dDaFile(LuRInv(2),iOpt,Wrk(kRia2+iOffset2),lTot,iAdr)
      end do
    end do

  end do    !iBat

  call ChoMP2_OpenF(iClos,iTypR,iSym)

end do !iSym

call DaClos(LuRInv(1))
call DaClos(LuRInv(2))

end subroutine ChoMP2g_Reord_R
