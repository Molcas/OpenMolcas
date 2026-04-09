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
! Copyright (C) 1993,1999, Roland Lindh                                *
!***********************************************************************

subroutine TnsCtl(Wrk,nWrk,nijkl,mabMax,mabMin,mcdMax,mcdMin,HMtrxAB,HMtrxCD,la,lb,lc,ld,iCmpa,jCmpb,kCmpc,lCmpd,iShlla,jShllb, &
                  kShllc,lShlld,i_out)
!***********************************************************************
!                                                                      *
! Object: to transform the intermediate integral set directly to       *
!         the final integral set.                                      *
!         Note that the position in memory of the final set is not     *
!         fixed.                                                       *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             Modified by R.L Februrary, 1999.                         *
!***********************************************************************

use Basis_Info, only: Shells
use Breit, only: nComp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nWrk, nijkl, mabMax, mabMin, mcdMax, mcdMin, la, lb, lc, ld, iCmpa, jCmpb, kCmpc, lCmpd, iShlla, &
                                 jShllb, kShllc, lShlld
real(kind=wp), intent(inout) :: Wrk(nWrk)
real(kind=wp), intent(in) :: HMtrxAB(*), HMtrxCD(*)
integer(kind=iwp), intent(out) :: i_out
integer(kind=iwp) :: i_In, iW2, iW3, nab, ncd, nDim, ne, nf, nfijkl, nijklab

! If nComp==1
!   Integral are stored as e,f,IJKL in Wrk
! If nComp/=1
!   Integral are stored as ncomp,e,f,IJKL in Wrk

! Observe that Transf is false for s and p functions.

ne = (mabMax-mabMin+1)
nf = (mcdMax-mcdMin+1)
nab = iCmpa*jCmpb
ncd = kCmpc*lCmpd
nDim = max(ne*nf,nab*nf,nab*ncd)

iW2 = 1
iW3 = 1+nijkl*nDim

if (nComp /= 1) then
  Wrk(iW3:iW3+ne*nf*nijkl-1) = Wrk(iW2:iW2+ne*nf*nijkl-1)
  call DGetMO(Wrk(iW3),nComp,nComp,ne*nf*(nijkl/nComp),Wrk(iW2),ne*nf*(nijkl/nComp))
end if

! If (ss|ss) integral exit.

if (la+lb+lc+ld == 0) then
  i_out = 1
  return
end if

! Transpose if no transformation is needed.

if ((la*lb == 0) .and. (lc*ld == 0) .and. (.not. Shells(iShlla)%Transf) .and. (.not. Shells(jShllb)%Transf) .and. &
    (.not. Shells(kShllc)%Transf) .and. (.not. Shells(lShlld)%Transf)) then
  call DGeTMO(Wrk(iW2),ne*nf,ne*nf,nijkl,Wrk(iW3),nijkl)
  i_out = iW3
  return
end if

! Form matrix corresponding to the transfer equation, le,la,lb.
! The matrix transforms directly from e0 cartesians to real
! spherical harmonics.

if (la+lb == 0) then
  i_in = iW2
  i_out = iW3
else if ((la*lb == 0) .and. (.not. Shells(iShlla)%Transf) .and. (.not. Shells(jShllb)%Transf)) then
  call DGeTMO(Wrk(iW2),ne,ne,nf*nijkl,Wrk(iW3),nf*nijkl)
  i_in = iW3
  i_out = iW2
else

  ! Now transform directly (e,[f,IJKL]) to ([f,IJKL],AB)
  ! Int(lf*IJKL,lA*lB)=Int(le,lf*IJKL)*HMtrx(le,lA*lB)

  nfijkl = nf*nijkl
  call Sp_Mlt(Wrk(iW2),ne,Wrk(iW3),nfijkl,HMtrxAB,nab)
  i_in = iW3
  i_out = iW2
end if

! Form matrix corresponding to the transfer equation, lf,lC,lD

if (lc+ld == 0) then
  i_out = i_in
else if ((lc*ld == 0) .and. (.not. Shells(kShllc)%Transf) .and. (.not. Shells(lShlld)%Transf)) then
  call DGeTMO(Wrk(i_in),nf,nf,nijkl*nab,Wrk(i_out),nijkl*nab)
else

  ! Now transform directly (f,[IJKL,AB]) to ([IJKL,AB],CD)
  ! Int(IJKL*lA*lB,lC*lD)=Int(lf,IJKL*lA*lB)*HMtrx(lf,lC*lD)

  nijklAB = nijkl*nab
  call Sp_Mlt(Wrk(i_in),nf,Wrk(i_out),nijklAB,HMtrxCD,ncd)
end if

end subroutine TnsCtl
