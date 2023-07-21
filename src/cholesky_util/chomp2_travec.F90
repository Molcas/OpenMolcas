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

use ChoArr, only: iRS2F
use ChoSwp, only: IndRed

implicit real*8(a-h,o-z)
real*8 VecAO(*), VecMO(*), COcc(*), CVir(*)
real*8 Scr(lScr)
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2.fh"
character*13 SecNam
parameter(SecNam='ChoMP2_TraVec')
real*8 Fac(0:1)
data Fac/0.5d0,1.0d0/
! Statement function
MulD2h(i,j) = ieor(i-1,j-1)+1

if ((iLoc < 2) .or. (iLoc > 3)) then
  write(6,*) SecNam,': illegal iLoc = ',iLoc
  call ChoMP2_Quit(SecNam,'iLoc out of bounds!',' ')
end if
iSyScr = MulD2h(iSyCho,iSyCO)
if (lScr < nT1AOT(iSyScr)) then
  write(6,*) SecNam,': insufficient scratch space lScr = ',lScr
  write(6,*) SecNam,': needed                          = ',nT1AOT(iSyScr)
  call ChoMP2_Quit(SecNam,'Insufficient scratch space',' ')
else
  call FZero(Scr,nT1AOT(iSyScr))
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
        Go To 998
      end if
    end do
998 iSymBe = iSymAl
    iSymi = MulD2h(iSymBe,iSyCO)

    jAlpha = iAlpha-iBas(iSymAl)
    jBeta = iBeta-iBas(iSymBe)

    AOVal = Fac(min(abs(iAlpha-iBeta),1))*VecAO(iAlBe)
    kOffAl = iT1AOT(iSymi,iSymAl)+nOcc(iSymi)*(jAlpha-1)
    kOffBe = iT1AOT(iSymi,iSymBe)+nOcc(iSymi)*(jBeta-1)
    do i=1,nOcc(iSymi)
      Scr(kOffAl+i) = Scr(kOffAl+i)+AOVal*COcc(kOffBe+i)
      Scr(kOffBe+i) = Scr(kOffBe+i)+AOVal*COcc(kOffAl+i)
    end do

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
        Go To 999
      end if
    end do
999 iSymBe = MulD2h(iSymAl,iSyCho)

    jAlpha = iAlpha-iBas(iSymAl)
    jBeta = iBeta-iBas(iSymBe)

    AOVal = VecAO(iAlBe)

    iSymi = MulD2h(iSymBe,iSyCO)
    kOffAl = iT1AOT(iSymi,iSymAl)+nOcc(iSymi)*(jAlpha-1)
    kOffBe = iT1AOT(iSymi,iSymBe)+nOcc(iSymi)*(jBeta-1)
    do i=1,nOcc(iSymi)
      Scr(kOffAl+i) = Scr(kOffAl+i)+AOVal*COcc(kOffBe+i)
    end do

    iSymi = MulD2h(iSymAl,iSyCO)
    kOffAl = iT1AOT(iSymi,iSymAl)+nOcc(iSymi)*(jAlpha-1)
    kOffBe = iT1AOT(iSymi,iSymBe)+nOcc(iSymi)*(jBeta-1)
    do i=1,nOcc(iSymi)
      Scr(kOffBe+i) = Scr(kOffBe+i)+AOVal*COcc(kOffAl+i)
    end do

  end do

end if

! Second half-transformation step:
! VecMO(a,i) = sum_alpha CVir(alpha,a)*Scr(i,alpha)
! -------------------------------------------------

do iSymi=1,nSym

  iSyma = MulD2h(iSymi,iSyCho)
  iSymAl = MulD2h(iSyma,iSyCV)

  nTotAl = nBas(iSymAl)
  nTota = nVir(iSyma)
  nToti = nOcc(iSymi)

  if ((nToti > 0) .and. (nTota > 0) .and. (nTotAl > 0)) then
    kOff1 = iAOVir(iSymAl,iSyma)+1
    kOff2 = iT1AOT(iSymi,iSymAl)+1
    kOff3 = iT1am(iSyma,iSymi)+1
    call DGEMM_('T','T',nVir(iSyma),nOcc(iSymi),nBas(iSymAl),1.0d0,CVir(kOff1),nTotAl,Scr(kOff2),nToti,0.0d0,VecMO(kOff3),nTota)
  end if

end do

end subroutine ChoMP2_TraVec
