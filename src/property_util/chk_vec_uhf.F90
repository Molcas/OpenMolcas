!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************

subroutine Chk_vec_UHF(FName,Lu,isUHF)
! routine returns isUHF based on information in INPORB

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: FName
integer(kind=iwp), intent(inout) :: Lu
integer(kind=iwp), intent(out) :: isUHF
integer(kind=iwp) :: iVer, jVer
character(len=80) :: LINE
logical(kind=iwp) :: Exists
character(len=*), parameter :: Location = 'Chk_vec_UHF'
#include "inporbfmt.fh"

Line = 'not defined yet'

call OpnFl(FName,Lu,Exists)
if (.not. Exists) then
  write(u6,*) 'RdVec: File ',trim(FName),' not found!'
  call Abend()
end if
rewind(LU)
! Check version!
read(LU,'(A80)',end=999,ERR=999) Line
iVer = 0
do jVer=1,mxVer
  if (Magic(jVer) == Line(1:len(Magic(jVer)))) iVer = jVer
end do

if (iVer == 0) then
  call SysWarnMsg(Location,'INPORB file in old format',' ')
  call SysPutsEnd()
  isUHF = 0
  close(Lu)
  return
end if
50 continue
read(LU,'(A80)',end=999,ERR=999) Line
if (Line(1:5) /= '#INFO') goto 50
! Now Do the real job
read(Lu,'(a)',end=999,err=999) Line
read(Lu,*,end=999,err=999) isUHF
close(Lu)

return

999 continue
call SysWarnFileMsg(Location,FName,'Error during reading INPORB\n',Line)
call Abend()

end subroutine Chk_vec_UHF
