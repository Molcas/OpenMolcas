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

subroutine Untested(label)

implicit none
character(len=*), intent(in) :: label
#ifndef _DEVEL_
character(len=8) :: isdev
#endif

call WarningMessage(1,label//';This code is untested or experimental, and should be carefully verified.')
#ifndef _DEVEL_
isdev = ''
call getenvf('MOLCAS_ISDEV',isdev)
if (isdev == '') call abend()
#endif

end subroutine Untested
