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
subroutine SymAd1(lOper,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,AOInt,iBas,jBas,nIC,iIC,SOInt,nSOInt,nOp)
!***********************************************************************
!                                                                      *
! Object: to transform the one-electon matrix elements from AO basis   *
!         to SO basis.                                                 *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '91                                              *
!***********************************************************************

use Index_Functions, only: nTri3_Elem
use Basis_Info, only: Shells
use Symmetry_Info, only: iChBas, iChTbl, iOper, nIrrep, Prmt
use SOAO_Info, only: iAOtSO
use Real_Spherical, only: iSphCr
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: lOper, iAng, jAng, iCmp, jCmp, iShell, jShell, iShll, jShll, iAO, jAO, iBas, jBas, nIC, nSOInt, &
                                 nOp(2)
real(kind=wp), intent(in) :: AOInt(iBas*jBas,iCmp,jCmp,nIC)
integer(kind=iwp), intent(inout) :: iIC
real(kind=wp), intent(inout) :: SOInt(iBas*jBas,nSOInt)
integer(kind=iwp) :: i1, i2, iChBs, ii, iIrrep, j1, j12, j2, jChBs, jIC(0:7), jj, jMx, kIC, lSO
real(kind=wp) :: pae, pbr, xa, xb

#ifdef _DEBUGPRINT_
write(u6,*) ' lOper=',lOper
write(u6,*) ' nSOInt=',nSOInt
call RecPrt(' In SymAd1: AOInt',' ',AOInt,iBas*jBas,iCmp*jCmp*nIC)
call RecPrt(' In SymAd1: SOInt',' ',SOInt,iBas*jBas,nSOInt)
write(u6,*) ' iIC=',iIC
#endif
do iIrrep=0,nIrrep-1
  jIC(iIrrep) = -999999999
  if (.not. btest(lOper,iIrrep)) cycle
  jIC(iIrrep) = iIC
  iIC = iIC+1
end do

ii = nTri3_Elem(iAng)
jj = nTri3_Elem(jAng)

lSO = 0
do j1=0,nIrrep-1
  xa = real(iChTbl(j1,nOp(1)),kind=wp)
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle
    iChBs = iChBas(ii+i1)
    if (Shells(iShll)%Transf) iChBs = iChBas(iSphCr(ii+i1))
    pae = Prmt(iOper(nOp(1)),iChBs)

    do j2=0,nIrrep-1
      j12 = ieor(j1,j2)

      if (.not. btest(lOper,j12)) cycle
      kIC = jIC(j12)
      xb = real(iChTbl(j2,nOp(2)),kind=wp)
      jMx = jCmp
      if ((iShell == jShell) .and. (j1 == j2)) jMx = i1

      do i2=1,jMx
        if (iAOtSO(jAO+i2,j2) < 0) cycle
        lSO = lSO+1
        jChBs = iChBas(jj+i2)
        if (Shells(jShll)%Transf) jChBs = iChBas(iSphCr(jj+i2))
        pbr = Prmt(iOper(nOp(2)),jChBs)
        SOInt(:,lSO) = SOInt(:,lSO)+xa*pae*xb*pbr*AOInt(:,i1,i2,kIC)
      end do

    end do

  end do
end do
if (lSO /= nSOInt) then
  call WarningMessage(2,'Error in SymAd1, lSO /= nSOInt')
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In SymAd1: SOInt',' ',SOInt,iBas*jBas,nSOInt)
#endif

return

end subroutine SymAd1
