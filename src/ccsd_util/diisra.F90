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

subroutine diisra(wrk,wrksize,diispoint,num,v1,v2,v3,v4)
! this routine reads num vectors of amplitudes
! from prepaired diispoint(1) - diispoint(num) files
!
! diispoint - stack of lun numbers (I)
! num       - number of vectors to be read (1-4) (I)
! v1        - map type of vector 1 (I)
! v2        - map type of vector 1 (I)
! v3        - map type of vector 1 (I)
! v4        - map type of vector 1 (I)
! if there is less than 4 vectors required use any map type

use ccsd_global, only: Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, diispoint(4), num
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(inout) :: v1, v2, v3, v4
integer(kind=iwp) :: lun, rc

if (num == 1) then
  lun = diispoint(1)
  call getmediate(wrk,wrksize,lun,v1,rc)
else if (num == 2) then
  lun = diispoint(1)
  call getmediate(wrk,wrksize,lun,v1,rc)
  lun = diispoint(2)
  call getmediate(wrk,wrksize,lun,v2,rc)
else if (num == 3) then
  lun = diispoint(1)
  call getmediate(wrk,wrksize,lun,v1,rc)
  lun = diispoint(2)
  call getmediate(wrk,wrksize,lun,v2,rc)
  lun = diispoint(3)
  call getmediate(wrk,wrksize,lun,v3,rc)
else if (num == 4) then
  lun = diispoint(1)
  call getmediate(wrk,wrksize,lun,v1,rc)
  lun = diispoint(2)
  call getmediate(wrk,wrksize,lun,v2,rc)
  lun = diispoint(3)
  call getmediate(wrk,wrksize,lun,v3,rc)
  lun = diispoint(4)
  call getmediate(wrk,wrksize,lun,v4,rc)
end if

return

end subroutine diisra
