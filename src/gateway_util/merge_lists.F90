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

subroutine Merge_Lists(Mode,nAt)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
character, intent(in) :: Mode
integer(kind=iwp), intent(in) :: nAt
integer(kind=iwp) :: i_1, i_2, i_P, i_R, iOff_1, iOff_2, iOff_3, iOff_Iter, ipCx_1, ipCx_2, ipCx_P, ipCx_R, ipEner_1, ipEner_2, &
                     ipEner_P, ipEner_R, ipGx_1, ipGx_2, ipGx_P, ipGx_R, iter_1, iter_2, iter_3, iter_P, iter_R, n1, n2
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: rList(:,:)
integer(kind=iwp), allocatable :: iList(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Get the sizes of the arrays on the runfile. Same size for both files.

call qpg_iArray('Slapaf Info 1',Found,n1)
call qpg_dArray('Slapaf Info 2',Found,n2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over both files and pick up the fields and pick up pointers,
! and interations counts.

call mma_allocate(iList,n1,2,label='iList')
call mma_allocate(rList,n2,2,label='rList')
call NameRun('RUNREAC')
call Get_iArray('Slapaf Info 1',iList(1,1),n1)
call Get_dArray('Slapaf Info 2',rList(1,1),n2)
i_R = 1
iter_R = iList(2,1)
ipEner_R = 1+iList(5,1)
ipCx_R = 1+iList(6,1)
ipGx_R = 1+iList(7,1)
call NameRun('RUNPROD')
call Get_iArray('Slapaf Info 1',iList(1,2),n1)
call Get_dArray('Slapaf Info 2',rList(1,2),n2)
i_P = 2
iter_P = iList(2,2)
ipEner_P = 1+iList(5,2)
ipCx_P = 1+iList(6,2)
ipGx_P = 1+iList(7,2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Translate into generic variables, _1=from, _2=to.

if (Mode == 'R') then
  i_1 = i_P
  iter_1 = iter_P
  ipEner_1 = ipEner_P
  ipCx_1 = ipCx_P
  ipGx_1 = ipGx_P

  i_2 = i_R
  iter_2 = iter_R
  ipEner_2 = ipEner_R
  ipCx_2 = ipCx_R
  ipGx_2 = ipGx_R
else
  i_1 = i_R
  iter_1 = iter_R
  ipEner_1 = ipEner_R
  ipCx_1 = ipCx_R
  ipGx_1 = ipGx_R

  i_2 = i_P
  iter_2 = iter_P
  ipEner_2 = ipEner_P
  ipCx_2 = ipCx_P
  ipGx_2 = ipGx_P
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Start moving stuff around

! Update iteration counter

iter_3 = iter_2+1
iList(2,i_2) = iter_3
iOff_1 = (iter_1-1)*3*nAt
iOff_2 = (iter_2-1)*3*nAt
iOff_3 = (iter_3-1)*3*nAt

! Move the last item(s) in the "to" file up one step.

rList(ipEner_2+iter_3-1,i_2) = rList(ipEner_2+iter_2-1,i_2)
rlist(ipCx_2+iOff_3:ipCx_2+iOff_3+3*nAt-1,i_2) = rList(ipCx_2+iOff_2:ipCx_2+iOff_2+3*nAt-1,i_2)
rlist(ipGx_2+iOff_3:ipGx_2+iOff_3+3*nAt-1,i_2) = rList(ipGx_2+iOff_2:ipGx_2+iOff_2+3*nAt-1,i_2)

! Copy the last item(s) in the "from" file into the
! second last position of the "to" file.

rList(ipEner_2+iter_2-1,i_2) = rList(ipEner_1+iter_1-1,i_1)
rlist(ipCx_2+iOff_2:ipCx_2+iOff_2+3*nAt-1,i_2) = rList(ipCx_1+iOff_1:ipCx_1+iOff_1+3*nAt-1,i_1)
rlist(ipGx_2+iOff_2:ipGx_2+iOff_2+3*nAt-1,i_2) = rList(ipGx_1+iOff_1:ipGx_1+iOff_1+3*nAt-1,i_1)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out the stuff on the appropriate run file.

if (Mode == 'R') then
  call NameRun('RUNREAC')
else
  call NameRun('RUNPROD')
end if

call Put_iArray('Slapaf Info 1',iList(1,i_2),n1)
call Put_dArray('Slapaf Info 2',rList(1,i_2),n2)
call qpg_iScalar('iOff_Iter',Found)
if (Found) then
  call Get_iScalar('iOff_Iter',iOff_Iter)
  call Put_iScalar('iOff_Iter',iOff_Iter+1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate the memory allocations.

call mma_deallocate(rList)
call mma_deallocate(iList)
!                                                                      *
!***********************************************************************
!                                                                      *
! Open the default run file.

call NameRun('RUNFILE')
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Merge_Lists
