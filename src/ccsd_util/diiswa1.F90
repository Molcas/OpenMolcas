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

subroutine diiswa1(wrk,wrksize,diispoint)
! this routine does:
! a) upgrade diispoint
! b) write new amplitudes T21,T22,T23,T13,T14 into proper position
!
! diispoint - array of lun's where N-1, N-2 .. amplitudes are stored (I/O)

#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
integer diispoint(1:4)
! help variables
integer lun1, rc, p

!1 upgrade diispoint

lun1 = diispoint(cycext)

do p=cycext-1,1,-1
  diispoint(p+1) = diispoint(p)
end do

diispoint(1) = lun1

!2 save aplitudes

!2.1 rewind lun1 file
lun1 = diispoint(1)
call filemanager(2,lun1,rc)

!2.2 write T21
call wrtmediate(wrk,wrksize,lun1,mapdt21,mapit21,rc)

!2.3 write T22
call wrtmediate(wrk,wrksize,lun1,mapdt22,mapit22,rc)

!2.4 write T23
call wrtmediate(wrk,wrksize,lun1,mapdt23,mapit23,rc)

!2.5 write T13
call wrtmediate(wrk,wrksize,lun1,mapdt13,mapit13,rc)

!2.6 write T14
call wrtmediate(wrk,wrksize,lun1,mapdt14,mapit14,rc)

!2.1 rewind lun1 file
call filemanager(2,lun1,rc)

return

end subroutine diiswa1
