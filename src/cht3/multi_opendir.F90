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

subroutine multi_opendir(FNAM,iunit)
! Direct fortran I/O with irregular data records
!
! Assume here RECL in byte units  (Same assumption in t3smat.f).
!
! PV/LAOG, 22 may 2003.

use ChT3_global, only: nblock
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: FNAM
integer(kind=iwp), intent(in) :: iunit
integer(kind=iwp) :: iost
logical(kind=iwp) :: is_error

!open(unit=iunit,file=FNAM,access='direct',form='unformatted',status='unknown',recl=nblock*8)
call MOLCAS_Open_Ext2(iUnit,FNam,'direct','unformatted',iost,.true.,nblock*8,'unknown',is_error)
if ((iost > 0) .or. is_error) write(u6,*) 'Multi_OpenDir: Error opening file!'

return

end subroutine multi_opendir
