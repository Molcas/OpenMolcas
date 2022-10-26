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

subroutine cct3_getint(wrk,wrksize,i,symi,r,rc)
! this routine reads integrals R_i(a,bc) for given i in given symi
!
! i    - number of orbital (I)
! symi - irrep of i (I)
! r    - R (I/O)
! rc   - return (error) code (O)

use CCT3_global, only: daddr, Map_Type, noa, T3IntPos, t3nam
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, i, symi
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
type(Map_Type), intent(inout) :: r
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp) :: iadd, im, isym, length, lun, num, pos !, rc1

!1 some tests

if (i > noa(symi)) then
  ! RC=1 : i is higher than occipied in this irrep (Stup)
  rc = 1
  return
end if

if (i < 1) then
  ! RC=2 : i is less than 1 (Stup)
  rc = 2
  return
end if

!2 calc number for this orbital

iadd = 0
if (symi > 1) then
  do isym=1,symi-1
    iadd = iadd+noa(isym)
  end do
end if

num = iadd+i

!3 get R

lun = 1
daddr(lun) = T3IntPos(num)

call daname(lun,t3nam)

call idafile(lun,2,r%d,513*6,daddr(lun))
call idafile(lun,2,r%i,8*8*8,daddr(lun))

pos = r%pos0
length = 0
do im=1,r%d(0,5)
  r%d(im,1) = pos
  pos = pos+r%d(im,2)
  length = length+r%d(im,2)
  !write(u6,99) ' MAP',(r%d(im,k),k=1,6)
end do

if (length > 0) call ddafile(lun,2,wrk(r%pos0),length,daddr(lun))
!call cct3_getmediate(wrk,wrksize,lun,r,rc1)

call daclos(lun)

return

!99 format(a3,i8,2x,i8,4(2x,i2))

end subroutine cct3_getint
