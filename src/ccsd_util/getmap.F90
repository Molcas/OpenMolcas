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

subroutine getmap(lun,length,map,rc)
! this routine reads %d and %i of the given mediate
! from lun and reconstructs %d to actual positions %pos0
!
! lun    - Logical unit number of file, where mediate is stored (Input)
! length - overall length of mediate (Output)
! map    - map type corresponding to given mediate (Output)
! rc     - return (error) code (Output)
!
! N.B.
! all mediates are stored as follows
! 1 - Map_Type
! 2 - one record with complete mediate

use ccsd_global, only: daddr, iokey, Map_Type
use Definitions, only: iwp
#ifdef __INTEL_COMPILER
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: lun
integer(kind=iwp), intent(out) :: length, rc
type(Map_Type), intent(inout) :: map
integer(kind=iwp) :: im, pos

rc = 0

!1 read %d

if (iokey == 1) then
  ! Fortran IO
  read(lun) map%d(:,:),map%i(:,:,:)
# ifdef __INTEL_COMPILER
  ! workaround for a bug in some Intel versions
else if (iokey < 0) then
  write(u6,*) 'this should never happen'
# endif

else
  ! MOLCAS IO
  call idafile(lun,2,map%d,size(map%d),daddr(lun))
  call idafile(lun,2,map%i,size(map%i),daddr(lun))
end if

!2 change positions in %d to proper one and calculate overall length

pos = map%pos0
length = 0

do im=1,map%d(0,5)

  map%d(im,1) = pos
  pos = pos+map%d(im,2)
  length = length+map%d(im,2)

end do

return

end subroutine getmap
