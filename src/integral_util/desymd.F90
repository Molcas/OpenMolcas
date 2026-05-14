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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DesymD(lOper,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,DAO,iBas,jBas,DSO,nDSO,nOp,FactNd)
!***********************************************************************
!                                                                      *
! Object: desymmetrize the first order density matrix.                 *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use Index_Functions, only: nTri3_Elem
use Basis_Info, only: Shells
use Symmetry_Info, only: iChBas, iChTbl, iOper, nIrrep, Prmt
use SOAO_Info, only: iAOtSO
use Real_Spherical, only: iSphCr
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: lOper, iAng, jAng, iCmp, jCmp, iShell, jShell, iShll, jShll, iAO, jAO, iBas, jBas, nDSO, nOp(2)
real(kind=wp), intent(out) :: DAO(iBas*jBas,iCmp,jCmp)
real(kind=wp), intent(in) :: DSO(iBas*jBas,nDSO), FactNd
integer(kind=iwp) :: i1, i2, iChBs, ii, j1, j12, j2, jChBs, jj, jMx, lSO
real(kind=wp) :: Deg, FactNs, pa, Xa, Xb

#ifdef _DEBUGPRINT_
write(u6,*) ' lOper=',lOper
call RecPrt(' In DesymD: DSO',' ',DSO,iBas*jBas,nDSO)
#endif

DAO(:,:,:) = Zero
lSO = 1
ii = nTri3_Elem(iAng)
jj = nTri3_Elem(jAng)
do j1=0,nIrrep-1
  Xa = real(iChTbl(j1,nOp(1)),kind=wp)
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle
    iChBs = iChBas(ii+i1)
    if (Shells(iShll)%Transf) iChBs = iChBas(iSphCr(ii+i1))
    pa = Prmt(iOper(nOp(1)),iChBs)

    do j2=0,j1
      j12 = ieor(j1,j2)
      if (.not. btest(lOper,j12)) cycle
      Xb = real(iChTbl(j2,nOp(2)),kind=wp)
      jMx = jCmp
      if ((iShell == jShell) .and. (j1 == j2)) jMx = i1
      do i2=1,jMx
        if (iAOtSO(jAO+i2,j2) < 0) cycle
        jChBs = iChBas(jj+i2)
        if (Shells(jShll)%Transf) jChBs = iChBas(iSphCr(jj+i2))

        Deg = Two
        if ((j1 == j2) .and. (iShell == jShell) .and. (i1 == i2)) Deg = One

        ! Parity factor due to symmetry operations applied to
        ! angular part of the basis function.
        FactNs = pa*Prmt(iOper(nOp(2)),jChBs)
        DAO(:,i1,i2) = DAO(:,i1,i2)+Deg*Xa*Xb*FactNs*DSO(:,lSO)
        lSO = lSO+1
      end do
    end do

  end do
end do

if (FactNd /= One) DAO(:,:,:) = FactNd*DAO(:,:,:)
#ifdef _DEBUGPRINT_
call RecPrt(' In DesymD: DAO',' ',DAO,iBas*jBas,iCmp*jCmp)
#endif

end subroutine DesymD
