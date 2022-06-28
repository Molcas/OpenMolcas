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

subroutine Cp_Prop_Int(A_Int,nAInt,B_Int,nBInt,nBas,nIrrep,Label)
!***********************************************************************
!                                                                      *
! Object: replace the diagonal blocks of the property integrals.       *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAInt, nBInt, nIrrep, nBas(0:nIrrep-1), Label
real(kind=wp), intent(out) :: A_Int(nAInt)
real(kind=wp), intent(in) :: B_Int(nBInt)
integer(kind=iwp) :: iCmp, iExp, iIrrep, iLen, jIrrep, Len_

iCmp = 1
iExp = 1
do iIrrep=1,nIrrep
  do jIrrep=1,iIrrep
    if (.not. btest(Label,Mul(iIrrep,jIrrep)-1)) cycle
    if (iIrrep == jIrrep) then
      Len_ = nBas(iIrrep-1)*(nBas(iIrrep-1)+1)/2
      do iLen=0,Len_-1
        !write(u6,*) A_Int(iExp+iLen),B_Int(iCmp+iLen)
        A_Int(iExp+iLen) = B_Int(iCmp+iLen)
      end do
      iCmp = iCmp+Len_
      iExp = iExp+Len_
    else
      Len_ = nBas(iIrrep-1)*nBas(jIrrep-1)
      iExp = iExp+Len_
    end if
  end do
end do

return

end subroutine Cp_Prop_Int
