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

subroutine cct3_getmap(lun,poss0,length,mapd,mapi,rc)
! this routine reads mapd and mapi of the given mediate
! from lun and reconstructs mapd to actual positions poss0
!
! lun    - Logical unit number of file, where mediate is stored (Input)
! poss0  - initial position in WRK, where mediate will be stored (Input)
! length - overall length of mediate (Output)
! mapd   - direct map matrix corresponding to given mediate (Output)
! mapi   - inverse map matrix corresponding to given mediate (Output)
! rc     - return (error) code (Output)
!
! N.B.
! all mediates are stored as follows
! 1 - mapd, mapi
! 2 - one record with complete mediate

#include "t31.fh"

integer lun, poss0, rc
integer mapd(0:512,1:6)
integer mapi(1:8,1:8,1:8)
! help variables
integer poss, im, length

rc = 0

!1 read mapd

if (iokey == 1) then
  ! Fortran IO
  read(lun) mapd,mapi

else
  ! MOLCAS IO
  call idafile(lun,2,mapd,513*6,daddr(lun))
  call idafile(lun,2,mapi,8*8*8,daddr(lun))
end if

!2 change positions in mapd to proper one and calculate overall length

poss = poss0
length = 0

do im=1,mapd(0,5)

  mapd(im,1) = poss
  poss = poss+mapd(im,2)
  length = length+mapd(im,2)

end do

return

end subroutine cct3_getmap
