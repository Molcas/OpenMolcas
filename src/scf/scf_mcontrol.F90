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

subroutine Scf_Mcontrol(id_call)

use Para_Info, only: MyRank
use InfSCF, only: nIter, DThr, EThr, FThr
use ChoSCF, only: ALGO, dmpk, nScreen
use Cholesky, only: timings
use Constants, only: Zero

implicit none
integer id_call
integer icount, icount0
character(len=512) List
character(len=32) value

icount = 0

if (id_call == 1) then

  ! Label definitions

  write(List,100) 'SCF_started_OK:(-:-):',ALGO,timings,dmpK,EThr,DThr,FThr,nIter(1),nScreen

  ! Initialize the control system

  call MolcasControlInit(List)
  return

else

  ! Read the molcas control file

  !1
  call MolcasControl('Cho_ALGO',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) ALGO
    write(6,*) '--- Warning: Cho_ALGOrithm changed by user to the value ',ALGO
    icount = icount+1
  end if
  !2
  call MolcasControl('Chotime',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) timings
    write(6,*) '--- Warning: Cholesky timings visualization changed by user to the value ',timings
    icount = icount+1
  end if
  !3
  call MolcasControl('En_thr',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) Ethr
    write(6,*) '--- Warning: SCF Energy threshold changed by user to the value ',Ethr
    icount = icount+1
  end if
  !4
  call MolcasControl('D_thr',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) Dthr
    write(6,*) '--- Warning: SCF Density threshold changed by user to the value ',Dthr
    icount = icount+1
  end if
  !5
  call MolcasControl('F_thr',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) Fthr
    write(6,*) '--- Warning: SCF Fmat threshold changed by user to the value ',Fthr
    icount = icount+1
  end if
  !6
  call MolcasControl('MaxIter',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) nIter(1)
    write(6,*) '--- Warning: SCF Max # iterations changed by user to the value ',nIter(1)
    icount = icount+1
  end if
  !7
  call MolcasControl('nScreen',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) nScreen
    write(6,*) '--- Warning: Cholesky LK option nSCREEN changed by user to the value ',nScreen
    icount = icount+1
  end if
  !8
  call MolcasControl('dmpK',value)
  if (value(1:4) /= '    ') then
    read(value,*,err=101,end=102) dmpK
    write(6,*) '--- Warning: Cholesky LK option DMPK changed by user to the value ',dmpK
    icount = icount+1
  end if

end if

icount0 = icount

! Get the true updated counter in parallel runs

call gaIgOP_SCAL(icount,'max')

if (MyRank == 0) then
  if (icount > icount0) then
    write(6,*) ' Steering will NOT be activated this time because'
    write(6,*) ' molcas.control file must be changed on node_0 !!'
    call gaIgOP_SCAL(icount,'min')
  end if
end if

if (icount > 0) then

  ! Trick to broadcast the values across nodes

  if (MyRank /= 0) then
    ALGO = 0
    nIter(1) = 0
    nScreen = 0
    dmpK = Zero
    EThr = Zero
    DThr = Zero
    FThr = Zero
  end if
  call gaIgOP_SCAL(ALGO,'+')
  call gaIgOP_SCAL(nScreen,'+')
  call gaIgOP(nIter(1),1,'+')
  call gadgOP_SCAL(dmpK,'+')
  call gadgOP_SCAL(EThr,'+')
  call gadgOP_SCAL(DThr,'+')
  call gadgOP_SCAL(FThr,'+')

  ! Update label values (note that "timings" is locally updated!)

  write(List,100) 'SCF_modified_by_user:',ALGO,timings,dmpK,EThr,DThr,FThr,nIter(1),nScreen

  ! Initialize the control system with the new values

  call MolcasControlInit(List)
  return

end if

return

100 format(A21,',Cho_ALGO=',I2,',Chotime=',L2,',dmpK=',ES11.4,',En_thr=',ES11.4,',D_thr=',ES11.4,',F_thr=',ES11.4,',MaxIter=',I4, &
           ',nScreen=',I4)
101 write(6,*) 'Scf_Mcontrol: error in data Input. ( icount= ',icount,' )'
102 write(6,*) 'Scf_Mcontrol: reached end of file. ( icount= ',icount,' )'

end subroutine Scf_Mcontrol
