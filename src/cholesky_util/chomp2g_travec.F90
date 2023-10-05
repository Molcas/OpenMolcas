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

subroutine ChoMP2g_TraVec(VecAO,VecMO,COrb1,COrb2,Scr,lScr,iSyCho,iSyCO,iSyCV,iLoc,iMoType1,iMoType2)
!
! Jonas Bostrom, Feb 2010
!
! Purpose: compute pq-vector from reduced set AO vector.

use Symmetry_Info, only: Mul
use Cholesky, only: iBas, iiBstR, IndRed, iRS2F, nBas, nnBstR, nSym
use ChoMP2, only: iAoMo, iMoAo, iMoMo, nMo, nMoAo, nMoType
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: VecAO(*), COrb1(*), COrb2(*)
real(kind=wp), intent(_OUT_) :: VecMO(*)
integer(kind=iwp), intent(in) :: lScr, iSyCho, iSyCO, iSyCV, iLoc, iMoType1, iMoType2
real(kind=wp), intent(out) :: Scr(lScr)
integer(kind=iwp) :: iAlBe, iAlpha, iBeta, iSym, iSymAl, iSymBe, iSymP, iSymq, iSyScr, iVecType, jAlBe, jAlpha, jBeta, kOff1, &
                     kOff2, kOff3, kOffAl, kOffBe, nTotAl, nTotp, nTotq
real(kind=wp) :: AOVal
real(kind=wp), parameter :: Fac(0:1) = [Half,One]
character(len=*), parameter :: SecNam = 'ChoMP2_TraVec'

! Check what type of Cholesky vector to make (fro-occ, occ-occ.....)
iVecType = iMoType2+(iMoType1-1)*nMoType

if ((iLoc < 2) .or. (iLoc > 3)) then
  write(u6,*) SecNam,': illegal iLoc = ',iLoc
  call SysAbendMsg(SecNam,'iLoc out of bounds!',' ')
end if
iSyScr = Mul(iSyCho,iSyCO)
if (lScr < nMoAo(iSyScr,iMoType1)) then
  write(u6,*) SecNam,': insufficient scratch space lScr = ',lScr
  write(u6,*) SecNam,': needed                          = ',nMoAo(iSyScr,iMoType1)
  call SysAbendMsg(SecNam,'Insufficient scratch space',' ')
else
  Scr(1:nMoAo(iSyScr,iMoType1)) = Zero
end if

! First half-transformation step:
! Scr(i,alpha) = sum_beta VecAO(alpha,beta)*COrb1(i,beta)
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
    iSymP = Mul(iSymBe,iSyCO)

    jAlpha = iAlpha-iBas(iSymAl)
    jBeta = iBeta-iBas(iSymBe)

    AOVal = Fac(min(abs(iAlpha-iBeta),1))*VecAO(iAlBe)
    kOffAl = iMoAo(iSymP,iSymAl,iMoType1)+nMo(iSymP,iMoType1)*(jAlpha-1)
    kOffBe = iMoAo(iSymP,iSymBe,iMoType1)+nMo(iSymP,iMoType1)*(jBeta-1)
    Scr(kOffAl+1:kOffAl+nMo(iSymP,iMoType1)) = Scr(kOffAl+1:kOffAl+nMo(iSymP,iMoType1))+ &
                                               AOVal*COrb1(kOffBe+1:kOffBe+nMo(iSymP,iMoType1))
    Scr(kOffBe+1:kOffBe+nMo(iSymP,iMoType1)) = Scr(kOffBe+1:kOffBe+nMo(iSymP,iMoType1))+ &
                                               AOVal*COrb1(kOffAl+1:kOffAl+nMo(iSymP,iMoType1))

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

    iSymP = Mul(iSymBe,iSyCO)
    kOffAl = iMoAo(iSymP,iSymAl,iMoType1)+nMo(iSymP,iMoType1)*(jAlpha-1)
    kOffBe = iMoAo(iSymP,iSymBe,iMoType1)+nMo(iSymP,iMoType1)*(jBeta-1)
    Scr(kOffAl+1:kOffAl+nMo(iSymP,iMoType1)) = Scr(kOffAl+1:kOffAl+nMo(iSymP,iMoType1))+ &
                                               AOVal*COrb1(kOffBe+1:kOffBe+nMo(iSymP,iMoType1))

    iSymP = Mul(iSymAl,iSyCO)
    kOffAl = iMoAo(iSymP,iSymAl,iMoType1)+nMo(iSymP,iMoType1)*(jAlpha-1)
    kOffBe = iMoAo(iSymP,iSymBe,iMoType1)+nMo(iSymP,iMoType1)*(jBeta-1)
    Scr(kOffBe+1:kOffBe+nMo(iSymP,iMoType1)) = Scr(kOffBe+1:kOffBe+nMo(iSymP,iMoType1))+ &
                                               AOVal*COrb1(kOffAl+1:kOffAl+nMo(iSymP,iMoType1))

  end do

end if

! Second half-transformation step:
! VecMO(q,p) = sum_alpha COrb2(alpha,q)*Scr(p,alpha)
! -------------------------------------------------

do iSymp=1,nSym

  iSymq = Mul(iSymp,iSyCho)
  iSymAl = Mul(iSymq,iSyCV)

  nTotAl = nBas(iSymAl)
  nTotq = nMo(iSymq,iMoType2)
  nTotp = nMo(iSymp,iMoType1)

  if ((nTotp > 0) .and. (nTotq > 0) .and. (nTotAl > 0)) then
    kOff1 = iAoMo(iSymAl,iSymq,iMoType2)+1
    kOff2 = iMoAo(iSymp,iSymAl,iMoType1)+1
    kOff3 = iMoMo(iSymq,iSymp,iVecType)+1
    call DGEMM_('T','T',nMo(iSymq,iMoType2),nMo(iSymp,iMoType1),nBas(iSymAl),One,COrb2(kOff1),nTotAl,Scr(kOff2),nTotp,Zero, &
                VecMO(kOff3),nTotq)
  end if

end do

end subroutine ChoMP2g_TraVec
