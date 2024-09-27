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

subroutine VrfMtrx(Label,lOper,nComp,ip,Matrix)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden, January '91                  *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Basis_Info, only: nBas
use Gateway_global, only: PrPrt
use Symmetry_Info, only: Mul, nIrrep
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: nComp, lOper(nComp), ip(nComp)
real(kind=wp), intent(in) :: Matrix(*)
integer(kind=iwp) :: iComp, iIrrep, ip1, iSmLbl, jIrrep, n2
real(kind=wp) :: VrfSum
character(len=80) :: Line
real(kind=wp), external :: DDot_

call Untested('VrfMtrx')
do iComp=1,nComp
  VrfSum = Zero
  ip1 = ip(iComp)
  iSmLbl = lOper(iComp)
  if (Prprt) iSmLbl = merge(1,0,btest(iSmLbl,0))
  do iIrrep=0,nIrrep-1
    if (nBas(iIrrep) <= 0) cycle
    do jIrrep=0,iIrrep
      if (nBas(jIrrep) <= 0) cycle
      if (.not. btest(iSmLbl,Mul(iIrrep+1,jIrrep+1)-1)) cycle
      if (iIrrep == jIrrep) then
        n2 = nTri_Elem(nBas(iIrrep))
        VrfSum = VrfSum+DDot_(n2,Matrix(ip1),1,Matrix(ip1),1)
        ip1 = ip1+n2
      else
        n2 = nBas(iIrrep)*nBas(jIrrep)
        VrfSum = VrfSum+DDot_(n2,Matrix(ip1),1,Matrix(ip1),1)
        ip1 = ip1+n2
      end if
    end do
  end do
  ! Add the nuclear contribution and the operator position.
  n2 = 4
  VrfSum = VrfSum+DDot_(n2,Matrix(ip1),1,Matrix(ip1),1)
  write(Line,'(A,I5)') Label,iComp
  call Add_info(Line,[VrfSum],1,8)
end do

end subroutine VrfMtrx
