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

subroutine Put_nadc(colgradmode,Grad,nGrad)

use RunFile_data, only: lw
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: colgradmode, nGrad
real(kind=wp), intent(in) :: Grad(nGrad)
integer(kind=iwp) :: iGo
character(len=lw) :: Label

select case (colgradmode)
  case (0)
    Label = 'GRAD'
  case (1)
    Label = 'Grad State1'
  case (2)
    Label = 'Grad State2'
  case (3)
    Label = 'NADC'
  case default
    write(u6,*) 'put_nadc: invalid colgradmode',colgradmode
    call Abend()
end select
call Put_dArray(Label,Grad,nGrad)

call Get_iScalar('Grad ready',iGO)
iGO = ibset(iGO,1)
call Put_iScalar('Grad ready',iGO)

return

end subroutine Put_nadc
