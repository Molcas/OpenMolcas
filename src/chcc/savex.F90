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

subroutine SaveX(X,length,Lun,LunName,keyopen,keyclose)
! this routine does
! 1) keyopen = 1 - open LunName file with Lun
!              2 - rewind Lun file
!              3 - open LunName file with Lun with ACCESS='append'
!           else - nothing (i.e) file is opened
! 2) write X of dimension length
! 3) keyclose= 1 - close Lun file
!           else - nothing

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: length, Lun, keyopen, keyclose
real(kind=wp), intent(in) :: X(length)
character(len=6) :: LunName

!1
if (keyopen == 1) then
  !open(unit=Lun,file=LunName,form='unformatted')
  call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
else if (keyopen == 2) then
  rewind(Lun)
else if (keyopen == 3) then
  !mp open(unit=Lun,file=LunName,form='unformatted',ACCESS='append')

  call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
  call append_file_u(Lun)

end if

!2
call wri_chcc(Lun,length,X)

!3
if (keyclose == 1) close(Lun)

return

end subroutine SaveX
