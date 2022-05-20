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

subroutine read_abdata()

use Definitions, only: iwp, u6

implicit none
#include "abtab.fh"
integer(kind=iwp) :: i, ipos, k, lu_abdata, nerr
character(len=8) :: key
logical(kind=iwp) :: found_abdata
character(len=*), parameter :: ABDATA_NAME = 'ABDATA'
integer(kind=iwp), external :: isFreeUnit

call f_Inquire(ABDATA_NAME,found_abdata)
if (.not. found_abdata) then
  call warningmessage(2,' the abdata file does not exist.')
  call abend()
end if
lu_abdata = isFreeUnit(22)
call molcas_open(lu_abdata,ABDATA_NAME)

do
  read(lu_abdata,'(a8)') key
  if (key == 'NTAB1, N') exit
end do
read(lu_abdata,*) ntab1,ntab2,maxdeg
nerr = 0
if (ntab2-ntab1+1 > mxsiz2) then
  call warningmessage(2,' mxsiz2 is too small in readab.')
  write(u6,*) ' recompile. needs mxsiz2=',ntab2-ntab1+1
  nerr = 1
end if
if (maxdeg > mxsiz1) then
  call warningmessage(2,' mxsiz1 is too small in readab.')
  write(u6,*) ' recompile. needs mxsiz1=',maxdeg
  nerr = 1
end if
if (nerr == 1) call abend()
ipos = 0
do i=ntab1,ntab2
  do
    read(lu_abdata,'(a8)') key
    if (key == 'TAB POIN') exit
  end do
  ipos = ipos+1
  read(lu_abdata,*) k,tvalue(ipos),p0(ipos)
  read(lu_abdata,*)
  read(lu_abdata,*) (atab(k,ipos),k=0,maxdeg)
  read(lu_abdata,*)
  read(lu_abdata,*) (btab(k,ipos),k=0,maxdeg)
end do

close(lu_abdata)

return

end subroutine read_abdata
