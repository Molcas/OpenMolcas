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

subroutine RasScf_Mcontrol(id_call)

use Fock_util_global, only: ALGO, dmpk, Nscreen
use Cholesky, only: timings
use Para_Info, only: MyRank
use rasscf_global, only: MaxIt, Thre, ThrSX, ThrTE
use Constants, only: Zero
use Definitions, only: u6

implicit none
integer id_call
integer iCount, iCOunt0
character(len=512) List
character(len=32) value

icount = 0

if (id_call == 1) then

  ! Label definitions

  write(List,100) 'RASSCF_started_OK:(-:-):',ALGO,timings,dmpK,nScreen,MaxIt,ThrE,ThrSX,ThrTE

  ! Initialize the control system

  call MolcasControlInit(List)
  return

else

  ! Read the molcas control file
  !1
  call MolcasControl('Cho_ALGO',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) ALGO
    write(u6,*) '--- Warning: Cho_ALGO changed by user to the value ',ALGO
    icount = icount+1
  end if
  !2
  call MolcasControl('Chotime',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) timings
    write(u6,*) '--- Warning: Cholesky timings visualization changed by user to the value ',timings
    icount = icount+1
  end if
  !3
  call MolcasControl('nScreen',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) nScreen
    write(u6,*) '--- Warning: Cholesky LK option nSCREEN changed by user to the value ',nScreen
    icount = icount+1
  end if
  !4
  call MolcasControl('dmpK',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) dmpK
    write(u6,*) '--- Warning: Cholesky LK option DMPK changed by user to the value ',dmpK
    icount = icount+1
  end if
  !5
  call MolcasControl('MaxIter',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) MaxIt
    write(u6,*) '--- Warning: MaxIt changed by user to the value ',MaxIt
    icount = icount+1
  end if
  !6
  call MolcasControl('ThrE',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) ThrE
    write(u6,*) '--- Warning: ThrE changed by user to the value ',ThrE
    icount = icount+1
  end if
  !7
  call MolcasControl('ThrSX',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) ThrSX
    write(u6,*) '--- Warning: ThrSX changed by user to the value ',ThrSX
    icount = icount+1
  end if
  !8
  call MolcasControl('ThrTE',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) ThrTE
    write(u6,*) '--- Warning: ThrTE changed by user to the value ',ThrTE
    icount = icount+1
  end if

end if

icount0 = icount

! Get the true updated counter in parallel runs

call gaIgOP_SCAL(icount,'max')

if (MyRank == 0) then
  if (icount > icount0) then
    write(u6,*) ' Steering will NOT be activated this time because'
    write(u6,*) ' molcas.control file must be changed on node_0 !!'
    call gaIgOP_SCAL(icount,'min')
  end if
end if

if (icount > 0) then

  ! Trick to broadcast the values across nodes

  if (MyRank /= 0) then
    ALGO = 0
    MaxIt = 0
    nScreen = 0
    dmpK = Zero
    ThrE = Zero
    ThrSX = Zero
    ThrTE = Zero
  end if
  call gaIgOP_SCAL(ALGO,'+')
  call gaIgOP_SCAL(nScreen,'+')
  call gaIgOP_SCAL(MaxIt,'+')
  call gadgOP_SCAL(dmpK,'+')
  call gadgOP_SCAL(ThrE,'+')
  call gadgOP_SCAL(ThrSX,'+')
  call gadgOP_SCAL(ThrTE,'+')

  ! Update label values (note that "timings" is updated locally)

  write(List,100) 'RASSCF_modified_by_user:',ALGO,timings,dmpK,nScreen,MaxIt,ThrE,ThrSX,ThrTE

  ! Initialize the control system with the new values

  call MolcasControlInit(List)
  return

end if

return

100 format(A24,',Cho_ALGO=',I2,',Chotime=',L2,',dmpK=',ES11.4,',nScreen=',I4,',MaxIter=',I4,',ThrE=',ES11.4,',ThrSX=',ES11.4, &
           ',ThrTE=',ES11.4)

101 write(u6,*) 'RasScf_Mcontrol: error in data Input. ( icount= ',icount,' )'
102 write(u6,*) 'RasScf_Mcontrol: reached end of file. ( icount= ',icount,' )'

end subroutine RasScf_Mcontrol
