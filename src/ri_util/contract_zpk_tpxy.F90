!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine contract_Zpk_Tpxy(Zpk,nZpk,Txy,nTxy,Scrt,nScrt,Diag,nDiag,nnP,nBas_Aux,nAdens,nAvec,nAct,nIrrep)

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nZpk, nTxy, nScrt, nDiag, nIrrep, nnP(0:nIrrep-1), nBas_Aux(0:nIrrep-1), nAdens, nAvec, nAct(0:nIrrep-1)
real(kind=wp) :: Zpk(nZpk,nAVec), Txy(nTxy,nAdens), Scrt(nScrt), Diag(nDiag,nADens)
integer(kind=iwp) :: i, iDen, ip, iSym, j, jp, jSym, k, kp, kSym, l, nCumnnP, nCumnnP2, nCumnnP3

!***********************************************************************

do l=1,nAVec
  iDen = 1
  if (l > 1) iDen = 2

  nCumnnP = 0
  nCumnnP2 = 0
  nCumnnP3 = 0
  do iSym=0,nIrrep-1
    do i=1,nBas_Aux(iSym)
      ip = nCumnnP2+(i-1)*nnP(iSym)
      do j=1,nnP(iSym)
        Scrt(j) = Zero
        do k=1,nnP(iSym)
          kp = nCumnnP3+(k-1)*nnP(iSym)
          Scrt(j) = Scrt(j)+sign(One,Diag(nCumnnP+k,iDen))*Zpk(k+ip,l)*Txy(j+kp,iDen)
        end do
      end do
      do j=1,nnP(iSym)
        Zpk(j+ip,l) = Scrt(j)
      end do
    end do
    ! Now correct for the 2 factor
    do i=1,nBas_Aux(iSym)
      ip = nCumnnP2+(i-1)*nnP(iSym)
      do jSym=0,nIrrep-1
        kSym = ieor(jsym,isym)
        if (kSym <= jSym) then
          do j=1,nAct(jSym)
            if (kSym == jSym) then
              jp = ip+j*(j-1)/2
              do k=1,j-1
                Zpk(jp+k,l) = Zpk(jp+k,l)*Half
              end do
            else
              jp = ip+(j-1)*nAct(kSym)
              do k=1,nAct(kSym)
                Zpk(jp+k,l) = Zpk(jp+k,l)*Half
              end do
            end if
          end do
          if (kSym == jSym) then
            ip = ip+nAct(jSym)*(nAct(jSym)+1)/2
          else
            ip = ip+nAct(jSym)*nAct(kSym)
          end if
        end if
      end do
    end do

    !call RecPrt('Zpk',' ',Zpk(nCumnnP2+1,l),nnP(iSym),nBas_Aux(iSym))
    nCumnnP = nCumnnP+nnP(iSym)
    nCumnnP2 = nCumnnP2+nnP(iSym)*nBas_Aux(iSym)
    nCumnnP3 = nCumnnP3+nnP(iSym)**2
  end do
end do

end subroutine contract_Zpk_Tpxy
