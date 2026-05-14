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

subroutine saverest2(lunrst,energy,niter,iokey,daddr)
! this routine saves restart informations:
! energy, niter
! to prepaired position in lunrst

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lunrst, niter, iokey
real(kind=wp), intent(in) :: energy
integer(kind=iwp), intent(inout) :: daddr
integer(kind=iwp) :: idum(1)
real(kind=wp) :: dum(1)

!1 write energy,niter
if (iokey == 1) then
  ! Fortran IO
  write(lunrst) energy,niter

else
  ! MOLCAS IO
  dum(1) = energy
  call ddafile(lunrst,1,dum,1,daddr)
  idum(1) = niter
  call idafile(lunrst,1,idum,1,daddr)
end if

return

end subroutine saverest2
