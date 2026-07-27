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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine OLagFro2(NBSQT,DPT2,FPT2,ERI,Scr)

use caspt2_module, only: NBAS, NDEL, NFRO, NISH, NSYM
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NBSQT
real(kind=wp), intent(in) :: DPT2(NBSQT)
real(kind=wp), intent(inout) :: FPT2(NBSQT)
real(kind=wp), intent(out) :: ERI(NBSQT), Scr(NBSQT)
integer(kind=iwp) :: iMO, iOrb, iSym, iSymA, iSymB, iSymI, iSymJ, jOrb, nFroI, nIshI, nOrbI
real(kind=wp) :: Scal, Val

!write(u6,*) 'FPT2 before frozen orbital'
!call sqprt(fpt2,nbast)
iMO = 1
isymi = 1
isymj = 1
isyma = 1
isymb = 1
do iSym=1,nSym
  nOrbI = nBas(iSym)-nDel(iSym)
  nFroI = nFro(iSym)
  nIshI = nIsh(iSym)
  !! Fpq = ((pq|rs)-1/2(pr|qs))*Drs
  do iOrb=1,nFroI
    do jOrb=nFroI+1,nFroI+nIshI
      Scal = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
      call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,ERI,Scr)
      FPT2(iMO:iMO+nOrbI*nOrbI-1) = FPT2(iMO:iMO+nOrbI*nOrbI-1)+Scal*ERI(1:nOrbI*nOrbI)
      call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,ERI,Scr)
      FPT2(iMO:iMO+nOrbI*nOrbI-1) = FPT2(iMO:iMO+nOrbI*nOrbI-1)-Half*Scal*ERI(1:nOrbI*nOrbI)
    end do
  end do

  !! Symmetrize FPT2
  do iOrb=1,nOrbI
    do jOrb=1,iOrb-1
      Val = (FPT2(iMO+iOrb-1+nOrbI*(jOrb-1))+FPT2(iMO+jOrb-1+nOrbI*(iOrb-1)))*Half
      FPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
      FPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
    end do
  end do
  iMO = iMO+nOrbI*nOrbI
end do
!write(u6,*) 'FPT2 after frozen orbital'
!call sqprt(fpt2,nbast)

end subroutine OLagFro2
