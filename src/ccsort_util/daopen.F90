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

subroutine daopen(fname,lun,reclen)
! this routine opens direct access file
!
! fname - name of the file A8 (I)
! lun   - logical unit number (I)
! reclen- record length in R8 (I)

use Definitions, only: iwp

implicit none
character(len=8), intent(in) :: fname
integer(kind=iwp), intent(in) :: lun, reclen
integer(kind=iwp) :: f_iostat, recln
logical(kind=iwp) :: is_error

#ifdef _DECAXP_
recln = reclen*2
#else
recln = reclen*8
#endif

call molcas_open_ext2(lun,fname,'direct','unformatted',f_iostat,.true.,recln,'unknown',is_error)
!open(unit=lun,file=fname,form='unformatted',access='direct',recl=recln)

return

end subroutine daopen
