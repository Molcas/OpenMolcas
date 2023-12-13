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

subroutine LU2DESC(Lu,Desc)

use Fast_IO, only: CtlBlk, MxFile, pDesc, pHndle
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: Lu
integer(kind=iwp), intent(inout) :: Desc
integer(kind=iwp) :: handle, n, nFile
integer(kind=iwp), external :: lu2handle

handle = lu2handle(Lu)

n = 1
do
  if (CtlBlk(pHndle,n) == handle) exit
  n = n+1
  if (n > MxFile) then
    return
  end if
end do
nFile = n
Desc = CtlBlk(pDesc,nFile)

return

end subroutine LU2DESC
