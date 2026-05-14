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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Energy_Prt(Caller,Job,iBatch)
!
! Thomas Bondo Pedersen, March 2005.
!
! Purpose: print progress reports for the MP2 energy evaluation.

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Caller
integer(kind=iwp), intent(in) :: Job, iBatch
real(kind=wp) :: CME_Time(2,2) = Zero, CPU, Ratio, Wall
character(len=*), parameter :: SecNam = 'ChoMP2_Energy_Prt', ThisNam = 'Energy_Prt'

if (Job == 0) then

  CME_Time(:,:) = Zero

  write(u6,'(/,4X,A,/,4X,A)') 'Evaluation of MP2 energy correction','==================================='
  write(u6,'(4X,A,A)') 'Evaluator: ',Caller

  write(u6,'(/,4X,A,/,4X,A,/,4X,A)') 'Batch      CPU       Wall    Ratio',' No.     seconds    seconds', &
                                     '----------------------------------'

  call xFlush(u6)

else if (Job == 1) then

  call CWTime(CME_Time(1,1),CME_Time(2,1))

  call xFlush(u6)

else if (Job == 2) then

  call CWTime(CME_Time(1,2),CME_Time(2,2))
  CPU = CME_Time(1,2)-CME_Time(1,1)
  Wall = CME_Time(2,2)-CME_Time(2,1)
  if (abs(Wall) < 1.0e-8_wp) then
    if (abs(CPU) < 1.0e-8_wp) then
      Ratio = One
    else
      Ratio = 1.0e15_wp
    end if
  else
    Ratio = CPU/Wall
  end if
  write(u6,'(I9,2(1X,F10.2),1X,F6.3)') iBatch,CPU,Wall,Ratio

  call xFlush(u6)

else if (Job == 3) then

  write(u6,'(4X,A)') '----------------------------------'

  call xFlush(u6)

else

  call SysAbendMsg(SecNam,'Input parameter "Job" is out of range',' ')

end if

end subroutine ChoMP2_Energy_Prt
