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
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Sort_mck(A,B,iBas,jBas,kBas,lBas,iCmp,jCmp,kCmp,lCmp,iBasO,jBasO,kBasO,lBasO,iCmpO,jCmpO,kCmpO,lCmpO,nVec,nop,iAng,ishll)
!***********************************************************************
!                                                                      *
!     This subroutine is a stupid solution on a easy problem, but it   *
!     should work and it doesnt take to much CPU time.                 *
!     eaw                                                              *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri3_Elem
use Basis_Info, only: Shells
use Real_Spherical, only: iSphCr
use Symmetry_Info, only: iChBas, iOper, Prmt
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iBas, jBas, kBas, lBas, iCmp, jCmp, kCmp, lCmp, iBasO, jBasO, kBasO, lBasO, iCmpO, jCmpO, kCmpO, &
                                 lCmpO, nVec, nop(4), iAng(4), ishll(4)
real(kind=wp), intent(in) :: A(iBas*jBas*kBas*lBas,iCmp,jCmp,kCmp,lCmp,nVec)
real(kind=wp), intent(out) :: B(kBasO*kCmpO,lBasO,lCmpO,iBasO,iCmpO,jBasO,jCmpO*nvec)
integer(kind=iwp) :: iC, ichbs, ii, ijkl, iVec, jB, jC, jChBs, jj, kB, kC, kChBs, kk, lB, lC, lChBs, ll
real(kind=wp) :: PrA, pRb, pTc, pTSd, qFctr, rp

ii = nTri3_Elem(iAng(1))
jj = nTri3_Elem(iAng(2))
kk = nTri3_Elem(iAng(3))
ll = nTri3_Elem(iAng(4))

rp = One
do iVec=1,nVec
  do iC=1,iCmp
    ichbs = ichbas(ii+ic)
    if (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+ic))
    PrA = Prmt(iOper(nOp(1)),iChBs)
    do jC=1,jCmp
      jChBs = iChBas(jj+jc)
      if (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+jc))
      pRb = Prmt(iOper(nOp(2)),jChBs)
      do kC=1,kCmp
        kChBs = iChBas(kk+kc)
        if (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+kc))
        pTc = Prmt(iOper(nOp(3)),kChBs)
        do lC=1,lCmp
          lChBs = iChBas(ll+lC)
          if (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+lc))
          pTSd = Prmt(iOper(nOp(4)),lChBs)
          qFctr = pTSd*pTc*pRb*pRa*rp
          !*EAW 970930

          ! Some machines don't support more than 7 indices, need to fool around to make it work

          !do iB=1,iBas
          !  do jB=1,jBas
          !    do lB=1,lBas
          !      B(1:kBas,kC,lB,lC,iB,iC,jB,jC,iVec) = qfctr*A(iB,jB,1:kBas,lB,iC,jC,kC,lC,iVec)
          !    end do
          !  end do
          !end do
          ijkl = 0
          do lB=1,lBas
            do kB=1,kBas
              do jB=1,jBas
                B(kB+(kC-1)*kBasO,lB,lC,1:iBas,iC,jB,jC+(iVec-1)*jCmpO) = qfctr*A(ijkl+1:ijkl+iBas,iC,jC,kC,lC,iVec)
                ijkl = ijkl+iBas
              end do
            end do
          end do

          !EAW 970930
        end do
      end do
    end do
  end do
end do

return

end subroutine Sort_mck
