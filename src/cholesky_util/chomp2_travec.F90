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

subroutine ChoMP2_TraVec(VecAO,VecMO,COcc,CVir,Scr,lScr,iSyCho,iSyCO,iSyCV,iLoc)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: compute ai-vector from reduced set AO vector.

use Symmetry_Info, only: Mul
use Cholesky, only: iBas, iiBstR, IndRed, iRS2F, nBas, nnBstR, nSym
use ChoMp2, only: iAOVir, iT1am, iT1AOT, nOcc, nT1AOT, nVir
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: VecAO(*), COcc(*), CVir(*)
real(kind=wp), intent(_OUT_) :: VecMO(*)
integer(kind=iwp), intent(in) :: lScr, iSyCho, iSyCO, iSyCV, iLoc
real(kind=wp), intent(out) :: Scr(lScr)
integer(kind=iwp) :: iAlBe, iAlpha, iBeta, iSym, iSyma, iSymAl, iSymBe, iSymi, iSyScr, jAlBe, jAlpha, jBeta, kOff1, kOff2, kOff3, &
                     kOffAl, kOffBe, nTota, nTotAl, nToti
real(kind=wp) :: AOVal
real(kind=wp), parameter :: Fac(0:1) = [Half,One]
character(len=*), parameter :: SecNam = 'ChoMP2_TraVec'

if ((iLoc < 2) .or. (iLoc > 3)) then
  write(u6,*) SecNam,': illegal iLoc = ',iLoc
  call SysAbendMsg(SecNam,'iLoc out of bounds!',' ')
end if
iSyScr = Mul(iSyCho,iSyCO)
if (lScr < nT1AOT(iSyScr)) then
  write(u6,*) SecNam,': insufficient scratch space lScr = ',lScr
  write(u6,*) SecNam,': needed                          = ',nT1AOT(iSyScr)
  call SysAbendMsg(SecNam,'Insufficient scratch space',' ')
else
  Scr(1:nT1AOT(iSyScr)) = Zero
end if

! First half-transformation step:
! Scr(i,alpha) = sum_beta VecAO(alpha,beta)*COcc(i,beta)
! ------------------------------------------------------

if (iSyCho == 1) then

  do iAlBe=1,nnBstR(iSyCho,iLoc)

    jAlBe = IndRed(iiBstR(iSyCho,iLoc)+iAlBe,iLoc)
    iAlpha = iRS2F(1,jAlBe)
    iBeta = iRS2F(2,jAlBe)

    iSymAl = 1
    do iSym=nSym,2,-1
      if (iAlpha > iBas(iSym)) then
        iSymAl = iSym
        exit
      end if
    end do
    iSymBe = iSymAl
    iSymi = Mul(iSymBe,iSyCO)

    jAlpha = iAlpha-iBas(iSymAl)
    jBeta = iBeta-iBas(iSymBe)

    AOVal = Fac(min(abs(iAlpha-iBeta),1))*VecAO(iAlBe)
    kOffAl = iT1AOT(iSymi,iSymAl)+nOcc(iSymi)*(jAlpha-1)
    kOffBe = iT1AOT(iSymi,iSymBe)+nOcc(iSymi)*(jBeta-1)
    Scr(kOffAl+1:kOffAl+nOcc(iSymi)) = Scr(kOffAl+1:kOffAl+nOcc(iSymi))+AOVal*COcc(kOffBe+1:kOffBe+nOcc(iSymi))
    Scr(kOffBe+1:kOffBe+nOcc(iSymi)) = Scr(kOffBe+1:kOffBe+nOcc(iSymi))+AOVal*COcc(kOffAl+1:kOffAl+nOcc(iSymi))

  end do

else

  do iAlBe=1,nnBstR(iSyCho,iLoc)

    jAlBe = IndRed(iiBstR(iSyCho,iLoc)+iAlBe,iLoc)
    iAlpha = iRS2F(1,jAlBe)
    iBeta = iRS2F(2,jAlBe)

    iSymAl = 1
    do iSym=nSym,2,-1
      if (iAlpha > iBas(iSym)) then
        iSymAl = iSym
        exit
      end if
    end do
    iSymBe = Mul(iSymAl,iSyCho)

    jAlpha = iAlpha-iBas(iSymAl)
    jBeta = iBeta-iBas(iSymBe)

    AOVal = VecAO(iAlBe)

    iSymi = Mul(iSymBe,iSyCO)
    kOffAl = iT1AOT(iSymi,iSymAl)+nOcc(iSymi)*(jAlpha-1)
    kOffBe = iT1AOT(iSymi,iSymBe)+nOcc(iSymi)*(jBeta-1)
    Scr(kOffAl+1:kOffAl+nOcc(iSymi)) = Scr(kOffAl+1:kOffAl+nOcc(iSymi))+AOVal*COcc(kOffBe+1:kOffBe+nOcc(iSymi))

    iSymi = Mul(iSymAl,iSyCO)
    kOffAl = iT1AOT(iSymi,iSymAl)+nOcc(iSymi)*(jAlpha-1)
    kOffBe = iT1AOT(iSymi,iSymBe)+nOcc(iSymi)*(jBeta-1)
    Scr(kOffBe+1:kOffBe+nOcc(iSymi)) = Scr(kOffBe+1:kOffBe+nOcc(iSymi))+AOVal*COcc(kOffAl+1:kOffAl+nOcc(iSymi))

  end do

end if

! Second half-transformation step:
! VecMO(a,i) = sum_alpha CVir(alpha,a)*Scr(i,alpha)
! -------------------------------------------------

do iSymi=1,nSym

  iSyma = Mul(iSymi,iSyCho)
  iSymAl = Mul(iSyma,iSyCV)

  nTotAl = nBas(iSymAl)
  nTota = nVir(iSyma)
  nToti = nOcc(iSymi)

  if ((nToti > 0) .and. (nTota > 0) .and. (nTotAl > 0)) then
    kOff1 = iAOVir(iSymAl,iSyma)+1
    kOff2 = iT1AOT(iSymi,iSymAl)+1
    kOff3 = iT1am(iSyma,iSymi)+1
    call DGEMM_('T','T',nVir(iSyma),nOcc(iSymi),nBas(iSymAl),One,CVir(kOff1),nTotAl,Scr(kOff2),nToti,Zero,VecMO(kOff3),nTota)
  end if

end do

end subroutine ChoMP2_TraVec
