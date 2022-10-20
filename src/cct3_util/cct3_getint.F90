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

subroutine cct3_getint(wrk,wrksize,i,symi,posr0,mapdr,mapir,rc)
! this routine reads integrals R_i(a,bc) for given i in given symi
!
! i     - number of orbital (I)
! symi  - irrep of i (I)
! posr0 - initial position of R (I)
! mapdr - direct map of R (I)
! mapir - inverse map of R (I)
! rc    - return (error) code (O)

use CCT3_global, only: daddr, noa
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, i, symi, posr0, mapdr(0:512,6), mapir(8,8,8), rc
real(kind=wp) :: wrk(wrksize)
#include "t3int.fh"
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

call idafile(lun,2,mapdr,513*6,daddr(lun))
call idafile(lun,2,mapir,8*8*8,daddr(lun))

pos = posr0
length = 0
do im=1,mapdr(0,5)
  mapdr(im,1) = pos
  pos = pos+mapdr(im,2)
  length = length+mapdr(im,2)
  !write(u6,99) ' MAP',(mapdr(im,k),k=1,6)
  !99 format(a3,i8,2x,i8,4(2x,i2))
end do

if (length > 0) then
  call ddafile(lun,2,wrk(posr0),length,daddr(lun))
end if
!call cct3_getmediate(wrk,wrksize,lun,posr0,mapdr,mapir,rc1)

call daclos(lun)

return

end subroutine cct3_getint
