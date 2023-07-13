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
! Copyright (C) Giovanni Ghigo                                         *
!***********************************************************************

subroutine Rd_UDIC_RowH(nInter,nRowH,mRowH)
!***********************************************************************
!                                                                      *
! Object: Reading the Internal Coordinates required for Numerical      *
!         estimation of single rows and columns of Hessian             *
! Called from: PrePro when nRowH > 0                                   *
! Author: Giovanni Ghigo, University of Torino, Italy                  *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nInter, nRowH
integer(kind=iwp), intent(out) :: mRowH(nRowH)
integer(kind=iwp) :: iLines, iRowH, j, kLines, Lu, Lu_UDIC
character(len=120) :: Temp
character(len=16) :: filnam
character(len=8) :: cLbl
character(len=8), allocatable :: Labels(:)
integer(kind=iwp), external :: IsFreeUnit

Lu = u6
Lu_UDIC = IsFreeUnit(91)
filnam = 'UDIC'
call molcas_open(Lu_UDIC,filnam)
rewind(Lu_UDIC)
mRowH(:) = 0

! Find begining of definitions of internal coordinates

call mma_allocate(Labels,nInter,Label='Labels')

do
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:4) == 'VARY') exit
end do
do iLines=1,nInter
  do
    read(Lu_UDIC,'(A)') Temp
    call UpCase(Temp)
    if (Temp(1:3) /= 'FIX') exit
  end do
  cLbl = '        '
  do j=1,120
    if (Temp(j:j) == ' ') exit
    cLbl(j:j) = Temp(j:j)
  end do
  Labels(iLines) = cLbl
  do
    if (index(Temp,'&') == 0) exit
    read(Lu_UDIC,'(A)') Temp
  end do
end do

read(Lu_UDIC,'(A)') Temp ! Skip ROWH
outer: do iRowH=1,nRowH
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  cLbl(1:8) = Temp(1:8)
  do kLines=1,nInter
    if (cLbl == Labels(kLines)) then
      mRowH(iRowH) = kLines
      cycle outer
    end if
  end do
  call WarningMessage(2,'Error in rd_udic')
  write(Lu,*) '**********************************************'
  write(Lu,*) ' ERROR: Undefined internal ROWH coordinate in '
  write(Lu,*) ' ',Temp(1:60)
  write(Lu,*) '**********************************************'
  call Quit_OnUserError()
end do outer
close(Lu_UDIC)

call mma_deallocate(Labels)

return

end subroutine Rd_UDIC_RowH
