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

!----------------------------------------------------------------------*
! This subroutine reads from the formatted output of mpprop the        *
! coordinates of the expansion centers.                                *
!----------------------------------------------------------------------*
subroutine Get_Centers(nAt,xyz)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "files_qmstat.fh"
#include "warnings.h"
dimension xyz(MxAt,MxAt,3)
character*13 TheLine
logical Exist

! Open the file
Lu = 40
Lu = IsFreeUnit(40)
call Opnfl('MPPROP',Lu,Exist)
if (.not. Exist) then
  write(6,*)
  write(6,*) ' Can not locate output file from MpProp. '
  call Quit(_RC_IO_ERROR_READ_)
end if
rewind(Lu)

! Read until you get standard line
do
  read(Lu,'(A)') TheLine
  if (TheLine == '* All centers') exit
end do
read(Lu,*) i

! Read atom centers.
do i=1,nAt
  read(Lu,'(A)') TheLine
  read(Lu,*) (xyz(i,i,k),k=1,3)
  do j=1,10
    read(Lu,'(A)') TheLine
  end do
end do

! Read bond centers.
do i=2,nAt
  do j=1,i-1
    read(Lu,'(A)') TheLine
    read(Lu,*) (xyz(i,j,k),k=1,3)
    do jj=1,10
      read(Lu,'(A)') TheLine
    end do
  end do
end do

! Square xyz for later convenience
do i=2,nAt
  do j=1,i-1
    xyz(j,i,1) = xyz(i,j,1)
    xyz(j,i,2) = xyz(i,j,2)
    xyz(j,i,3) = xyz(i,j,3)
  end do
end do

close(Lu)

return

end subroutine Get_Centers
