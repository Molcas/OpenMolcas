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
! Copyright (C) 2001-2005, Valera Veryazov                             *
!***********************************************************************

subroutine molcas_binaryopen_vanilla(Lu,f_name)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: Lu
character(len=*), intent(in) :: f_name
integer(kind=iwp) :: lRealName
character(len=4096) :: RealName

!RealName = f_Name
call PrgmTranslate(f_Name,RealName,lRealName)
!write(u6,*) 'DEBUG binopen ',RealName(1:lRealName)
open(unit=Lu,file=RealName(1:lRealName),form='unformatted')

return

end subroutine molcas_binaryopen_vanilla
