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

subroutine Sort_mck(A,B,iBas,jBas,kBas,lBas,iCmp,jCmp,kCmp,lCmp,iBasO,jBasO,kBasO,lBasO,iCmpO,jCmpO,kCmpO,lCmpO,nVec,nop,iAng, &
                    indgrd,indgrd2,ishll,C)
!***********************************************************************
!                                                                      *
!     This subroutine is a stupid solution on a easy problem, but it   *
!     should work and it doesnt take to much CPU time.                 *
!     eaw                                                              *
!                                                                      *
!***********************************************************************

use Basis_Info, only: Shells
use Real_Spherical, only: iSphCr
use Symmetry_Info, only: iChBas, iOper, nIrrep
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iBas, jBas, kBas, lBas, iCmp, jCmp, kCmp, lCmp, iBasO, jBasO, kBasO, lBasO, iCmpO, jCmpO, kCmpO, lCmpO, nVec, &
                     nop(4), iAng(4), indgrd(3,4,0:nirrep-1), indgrd2(3,4,0:nirrep-1), ishll(4)
real(kind=wp) :: A(iBas*jBas*kBas*lBas,iCmp,jCmp,kCmp,lCmp,nVec), B(kBasO*kcmpO,lBasO,lcmpO,iBasO,iCmpO,jBasO,jCmpO*nvec), C(*)
integer(kind=iwp) :: iB, iC, ichbs, ii, ijkl, iVec, jB, jC, jChBs, jj, kB, kC, kChBs, kk, lB, lC, lChBs, ll
real(kind=wp) :: PrA, pRb, Prmt(0:7), pTc, pTSd, qFctr, rp
data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
! Statement Function
real(kind=wp) :: xPrmt
integer(kind=iwp) :: i, j, iOff, ixyz
xPrmt(i,j) = Prmt(iand(i,j))
iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6

ii = iOff(iAng(1))
jj = iOff(iAng(2))
kk = iOff(iAng(3))
ll = iOff(iAng(4))

rp = One
do iVec=1,nVec
  do iC=1,iCmp
    ichbs = ichbas(ii+ic)
    if (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+ic))
    PrA = xPrmt(iOper(nOp(1)),iChBs)
    do jC=1,jCmp
      jChBs = iChBas(jj+jc)
      if (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+jc))
      pRb = xPrmt(iOper(nOp(2)),jChBs)
      do kC=1,kCmp
        kChBs = iChBas(kk+kc)
        if (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+kc))
        pTc = xPrmt(iOper(nOp(3)),kChBs)
        do lC=1,lCmp
          lChBs = iChBas(ll+lC)
          if (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+lc))
          pTSd = xPrmt(iOper(nOp(4)),lChBs)
          qFctr = pTSd*pTc*pRb*pRa*rp
          !*EAW 970930

          ! Some machines don't support more than 7 indices, need to fool around to make it work

          !do iB=1,iBas
          !  do jB=1,jBas
          !    do kB=1,kBas
          !      do lB=1,lBas
          !        B(kB,kC,lB,lC,iB,iC,jB,jC,iVec) = qfctr*A(iB,jB,kB,lB,iC,jC,kC,lC,iVec)
          !      end do
          !    end do
          !  end do
          !end do
          ijkl = 0
          do lB=1,lBas
            do kB=1,kBas
              do jB=1,jBas
                do iB=1,iBas
                  ijkl = ijkl+1
                  B(kB+(kC-1)*kbaso,lB,lC,iB,iC,jB,jC+(iVec-1)*jcmpO) = qfctr*A(ijkl,iC,jC,kC,lC,iVec)
                end do
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
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(indgrd)
  call Unused_integer_array(indgrd2)
  call Unused_real_array(C)
end if

end subroutine Sort_mck
