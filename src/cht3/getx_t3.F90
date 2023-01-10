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

subroutine GetX_t3(X,length,Lun,LunName,keyopen,keyclose)
! this routine does:
! 1) keyopen = 0 - nothing (i.e) file is opened
!              1 - open LunName file with Lun
!              2 - rewind Lun file
!              3 - open LunName file with Lun with ACCESS='append'
! 2) read X of dimension length
! 3) keyclose = 0 - nothing
!               1 - close Lun file

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: length, Lun, keyopen, keyclose
real(kind=wp), intent(out) :: X(length)
character(len=6) :: LunName

!1
if (keyopen == 1) then
  !open(unit=Lun,file=LunName,form='unformatted')
  call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
else if (keyopen == 2) then
  rewind(Lun)
else if (keyopen == 3) then
  !mp !open(unit=Lun,file=LunName,form='unformatted',ACCESS='append')

  call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
  call append_file_u(Lun)

end if

!2
read(Lun) X(:)

!3
if (keyclose == 1) then
  close(Lun)
end if

return

end subroutine GetX_t3
