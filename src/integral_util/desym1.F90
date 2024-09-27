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
subroutine Desym1(lOper,iCmp,jCmp,iShell,jShell,iAO,jAO,DAO,iBas,jBas,DSO,nDSO,nOp,Scrt)
!***********************************************************************
!                                                                      *
! Object: desymmetrize the first order density matrix.                 *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChTbl
use SOAO_Info, only: iAOtSO
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: lOper, iCmp, jCmp, iShell, jShell, iAO, jAO, iBas, jBas, nDSO, nOp(2)
real(kind=wp), intent(out) :: DAO(iBas*jBas,iCmp,jCmp), Scrt(iBas*jBas)
real(kind=wp), intent(in) :: DSO(iBas*jBas,nDSO)
integer(kind=iwp) :: i1, i2, j1, j12, j2, jMx, lSO
real(kind=wp) :: Deg, Xa, Xb

#ifdef _DEBUGPRINT_
write(u6,*) ' lOper=',lOper
call RecPrt(' In Desym1: DSO',' ',DSO,iBas*jBas,nDSO)
#endif

DAO(:,:,:) = Zero

! D(P,Q)_ij = Sum(iSym,jSym) X(iSym,P) X(jSym,Q) D(iSym,jSym)_ij

! Loop over DSO, iSym >= jSym

lSO = 0
do j1=0,nIrrep-1
  Xa = real(iChTbl(j1,nOp(1)),kind=wp)
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle

    do j2=0,j1
      j12 = ieor(j1,j2)
      if (.not. btest(lOper,j12)) cycle
      Xb = real(iChTbl(j2,nOp(2)),kind=wp)
      jMx = jCmp
      if ((iShell == jShell) .and. (j1 == j2)) jMx = i1
      do i2=1,jMx
        if (iAOtSO(jAO+i2,j2) < 0) cycle
        lSO = lSO+1

        Deg = Two
        if (j1 == j2) Deg = One

        ! Parity factor due to symmetry operations applied to
        ! angular part of the basis function.

        DAO(:,i1,i2) = DAO(:,i1,i2)+Deg*Xa*Xb*DSO(:,lSO)

        if ((iShell == jShell) .and. (j1 == j2) .and. (i1 /= i2)) then
          call DGeTMO(DSO(1,lSO),iBas,iBas,jBas,Scrt,jBas)
          DAO(:,i2,i1) = DAO(:,i2,i1)+Deg*Xa*Xb*Scrt(:)
        end if
      end do
    end do

  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt(' In Desym1: DAO',' ',DAO,iBas*jBas,iCmp*jCmp)
#endif

return

end subroutine Desym1
