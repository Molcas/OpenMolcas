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

subroutine CmpInt(XInt,n_Int,nBas,nIrrep,Label)
!***********************************************************************
!                                                                      *
! Object: to remove the off-diagonal nonzero blocks of matrix elements *
!         for an operator.                                             *
!                                                                      *
!         XInt(n_Int):array with nonzero elements                      *
!                                                                      *
!         nBas(0:nIrrep-1):number of basis functions in each irrep     *
!                                                                      *
!         Label: symmetry label of the operator for which the          *
!                matrix elements where computed.                       *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             March 1991                                               *
!***********************************************************************

use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n_Int, nIrrep, nBas(0:nIrrep-1), Label
real(kind=wp) :: XInt(n_Int+4)
integer(kind=iwp) :: iadd, iCmp, iExp, iIrrep, ij, iLen, jIrrep, Len_

iCmp = 1
iExp = 1
do iIrrep=0,nIrrep-1
  do jIrrep=0,iIrrep
    ij = Mul(iIrrep+1,jIrrep+1)-1
    if (.not. btest(Label,ij)) cycle
    if (iIrrep == jIrrep) then
      Len_ = nBas(iIrrep)*(nBas(iIrrep)+1)/2
      do iLen=0,Len_-1
        XInt(iLen+iCmp) = Xint(iLen+iExp)
      end do
      iCmp = iCmp+Len_
      iExp = iExp+Len_
    else
      Len_ = nBas(iIrrep)*nBas(jIrrep)
      iExp = iExp+Len_
    end if
  end do
end do
do iadd=0,3
  XInt(iCmp+iadd) = XInt(iExp+iadd)
end do
n_Int = iCmp-1

return

end subroutine CmpInt
