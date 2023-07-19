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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine SlStatus(kIter,Energy,rGrad,Ex,nLines,delE,HUpMet,Step_Trunc,Print_Status)
!***********************************************************************
!                                                                      *
! Object:                                                              *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May '91                                                  *
!***********************************************************************

use Slapaf_Info, only: GrdLbl, GrdMax, iNeg, StpLbl, StpMax, UpMeth
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: kIter, nLines
real(kind=wp), intent(in) :: Energy, rGrad, Ex, delE
character(len=8), intent(in) :: HUpMet
character, intent(in) :: Step_Trunc
logical(kind=iwp), intent(in) :: Print_Status
#include "print.fh"
integer(kind=iwp) :: i, iPrint, iRout, iter, ivv, Lu, Lu_file, Lu_out, nvv
character(len=8) :: lNeg
character(len=128), allocatable :: Lines(:)
integer(kind=iwp), external :: isfreeunit

!                                                                      *
!***********************************************************************
!                                                                      *
Lu = u6
iRout = 52
iPrint = nPrint(iRout)

call mma_allocate(Lines,[-1,nLines],Label='Lines')

! Pick up previous energy

if (kIter == 1) then
  Lines(:) = ' '
  iter = 1
  write(Lines(-1),'(A)') '                       Energy     Grad      Grad              Step                 Estimated   '// &
                         'Geom       Hessian'
  write(Lines(0),'(A)') 'Iter      Energy       Change     Norm      Max    Element    Max     Element     Final Energy Update '// &
                        'Update   Index'
else
  call Get_cArray('Slapaf Info 3',Lines,(nLines+2)*128)
  ! Find first blank line
  iter = kIter
end if

if (iter > nLines) then
  call WarningMessage(2,'Status: iter > nLines')
  write(Lu,*) 'iter,nLines=',iter,nLines
  call abend()
else if (iter < 1) then
  call WarningMessage(2,'Status: iter < 1')
  call abend()
end if

write(lNeg,'(I3)') iNeg(1)
if (iNeg(2) /= iNeg(1)) then
  if (iNeg(2) > 99) then
    write(lNeg(4:8),'("(",I3,")")') iNeg(2)
  else if (iNeg(2) > 9) then
    write(lNeg(4:7),'("(",I2,")")') iNeg(2)
  else
    write(lNeg(4:6),'("(",I1,")")') iNeg(2)
  end if
end if
write(Lines(iter),100) iter,Energy,delE,rGrad,GrdMax,GrdLbl,StpMax,Step_Trunc,StpLbl,Ex,UpMeth,HUpMet,lNeg
!                                                                      *
!***********************************************************************
!                                                                      *
Lu_file = isfreeunit(8)

! Turn off updating the structure file during numerical integrations steps.

nvv = 1
if (Print_Status) nvv = 2
do ivv=1,nvv
  if (ivv == 1) then
    Lu_out = Lu
  else
    call Molcas_Open(Lu_File,'STRUCTURE')
    Lu_out = Lu_file
  end if

  if ((ivv == 2) .or. (iPrint >= 5)) then
    write(Lu_out,*)
    write(Lu_out,'(A)') '******************************************************************************************************'// &
                        '****************'
    write(Lu_out,'(A)') '*                                    Energy Statistics for Geometry Optimization                      '// &
                        '               *'
    write(Lu_out,'(A)') '******************************************************************************************************'// &
                        '****************'
    do i=-1,iter
      write(Lu_out,'(A118)') Lines(i)
    end do
    write(Lu_out,*)
  end if
  if (ivv == 2) close(Lu_file)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call Put_cArray('Slapaf Info 3',Lines(-1),(nLines+2)*128)
call mma_deallocate(Lines)

return

100 format(I3,F16.8,F12.8,2(F9.6,1X),A8,F9.6,A1,1X,A8,F16.8,1X,A,1X,A,1X,A)

end subroutine SlStatus
