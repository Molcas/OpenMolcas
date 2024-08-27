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

use Basis_Info, only: Shells
use Symmetry_Info, only: iChBas, iChTbl, iOper, nIrrep, Prmt
use SOAO_Info, only: iAOtSO
use Real_Spherical, only: iSphCr
use Constants, only: Zero, One, Two

implicit none
integer lOper, iAng, jAng, iCmp, jCmp, iShell, jShell, iShll, jShll, iAO, jAO, iBas, jBas, nDSO
real*8 DAO(iBas*jBas,iCmp,jCmp), DSO(iBas*jBas,nDSO)
integer nOp(2)
real*8 FactNd
integer lSO, ii, jj, j1, i1, j2, j12, i2, iChBs, jMx, jChBs
real*8 Xa, pa, Xb, FactNs, Deg

#ifdef _DEBUGPRINT_
write(6,*) ' lOper=',lOper
call RecPrt(' In DesymD: DSO',' ',DSO,iBas*jBas,nDSO)
#endif

call dcopy_(iBas*jBas*iCmp*jCmp,[Zero],0,DAO,1)
lSO = 1
ii = iAng*(iAng+1)*(iAng+2)/6
jj = jAng*(jAng+1)*(jAng+2)/6
do j1=0,nIrrep-1
  Xa = dble(iChTbl(j1,nOp(1)))
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle
    iChBs = iChBas(ii+i1)
    if (Shells(iShll)%Transf) iChBs = iChBas(iSphCr(ii+i1))
    pa = Prmt(iOper(nOp(1)),iChBs)

    do j2=0,j1
      j12 = ieor(j1,j2)
      if (iand(lOper,2**j12) == 0) Go To 300
      Xb = dble(iChTbl(j2,nOp(2)))
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
        call DaXpY_(iBas*jBas,Deg*Xa*Xb*FactNs,DSO(1,lSO),1,DAO(1,i1,i2),1)
        lSO = lSO+1
      end do
300   continue
    end do

  end do
end do

if (FactNd /= One) call DScal_(iBas*jBas*iCmp*jCmp,FactNd,DAO,1)
#ifdef _DEBUGPRINT_
call RecPrt(' In DesymD: DAO',' ',DAO,iBas*jBas,iCmp*jCmp)
#endif

end subroutine DesymD
