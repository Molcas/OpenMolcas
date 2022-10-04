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

subroutine Get_Grad(Grad,nGrad)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nGrad
real(kind=wp) :: Grad(nGrad)
integer(kind=iwp) :: mGrad
logical(kind=iwp) :: Found
character(len=16), parameter :: Label = 'GRAD'

call qpg_dArray(Label,Found,mGrad)
if ((.not. Found) .or. (nGrad == 0)) call SysAbendmsg('get_grad','Did not find:',Label)
if (mGrad /= nGrad) then
  write(u6,*) 'mGrad=',mGrad
  write(u6,*) 'nGrad=',nGrad
  call SysAbendmsg('get_grad','mGrad/=nGrad:',Label)
end if
call Get_dArray(Label,Grad,nGrad)

return

end subroutine Get_Grad
