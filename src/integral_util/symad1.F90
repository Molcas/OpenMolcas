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

use Basis_Info, only: Shells
use Symmetry_Info, only: iChBas, iChTbl, iOper, nIrrep, Prmt
use SOAO_Info, only: iAOtSO
use Real_Spherical, only: iSphCr

implicit none
integer lOper, iAng, jAng, iCmp, jCmp, iShell, jShell, iShll, jShll, iAO, jAO, iBas, jBas, nIC, iIC, nSOInt
integer nOp(2)
real*8 AOInt(iBas*jBas,iCmp,jCmp,nIC), SOInt(iBas*jBas,nSOInt)
integer jIC(0:7)
integer :: iTwoj(0:7) = [1,2,4,8,16,32,64,128]
integer iIrrep, ii, jj, lSO, j1, i1, iChBs, j2, j12, kIC, jMx, i2, jChBs
real*8 xa, pae, xb, pbr

#ifdef _DEBUGPRINT_
write(6,*) ' lOper=',lOper
write(6,*) ' nSOInt=',nSOInt
call RecPrt(' In SymAd1: AOInt',' ',AOInt,iBas*jBas,iCmp*jCmp*nIC)
call RecPrt(' In SymAd1: SOInt',' ',SOInt,iBas*jBas,nSOInt)
write(6,*) ' iIC=',iIC
#endif
do iIrrep=0,nIrrep-1
  jIC(iIrrep) = -999999999
  if (iand(lOper,iTwoj(iIrrep)) == 0) cycle
  jIC(iIrrep) = iIC
  iIC = iIC+1
end do

ii = iAng*(iAng+1)*(iAng+2)/6
jj = jAng*(jAng+1)*(jAng+2)/6

lSO = 0
do j1=0,nIrrep-1
  xa = dble(iChTbl(j1,nOp(1)))
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle
    iChBs = iChBas(ii+i1)
    if (Shells(iShll)%Transf) iChBs = iChBas(iSphCr(ii+i1))
    pae = Prmt(iOper(nOp(1)),iChBs)

    do j2=0,nIrrep-1
      j12 = ieor(j1,j2)

      if (iand(lOper,iTwoj(j12)) == 0) cycle
      kIC = jIC(j12)
      xb = dble(iChTbl(j2,nOp(2)))
      jMx = jCmp
      if ((iShell == jShell) .and. (j1 == j2)) jMx = i1

      do i2=1,jMx
        if (iAOtSO(jAO+i2,j2) < 0) cycle
        lSO = lSO+1
        jChBs = iChBas(jj+i2)
        if (Shells(jShll)%Transf) jChBs = iChBas(iSphCr(jj+i2))
        pbr = Prmt(iOper(nOp(2)),jChBs)
        call DaXpY_(iBas*jBas,xa*pae*xb*pbr,AOInt(1,i1,i2,kIC),1,SOInt(1,lSO),1)
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
