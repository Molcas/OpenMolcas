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

subroutine Start_Alaska()

use Slapaf_Info, only: iState, Request_Alaska, Request_RASSI
use UnixInfo, only: ProgName
use Definitions, only: iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: iEnd, istatus, LuInput, LuSpool, NACstatesOpt(2)
logical(kind=iwp) :: Exists
character(len=len(ProgName)) :: PName
character(len=180) :: Line
character(len=128) :: FileName
character(len=16) :: StdIn, JOB1, JOB2
character(len=8) :: Method
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
! Get the name of the module

PName = ProgName
call Upcase(PName)
PName = adjustl(PName)

iEnd = 1
do
  if (PName(iEnd:iEnd) == ' ') exit
  iEnd = iEnd+1
end do
iEnd = min(iEnd-1,5)
FileName = PName(1:iend)//'INP'

LuInput = 11
LuInput = IsFreeUnit(LuInput)
call StdIn_Name(StdIn)
call Molcas_Open(LuInput,StdIn)

if (Request_RASSI) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Request computation of overlaps.

  if (nPrint(1) >= 6) then
    write(u6,*)
    write(u6,*) ' Slapaf requests the computation of overlaps first!'
    write(u6,*)
  end if

  call Get_cArray('Relax Method',Method,8)
  if ((Method == 'CASPT2') .or. (Method == 'RASPT2')) then
    JOB1 = 'JOBMIX'
  else
    JOB1 = 'JOBIPH'
  end if
  call f_inquire('JOBAUTO',Exists)
  JOB2 = JOB1
  if (Exists) JOB2 = 'JOBAUTO'

  write(LuInput,'(A)') '>ECHO OFF'
  write(LuInput,'(A)') '> export SL_OLD_TRAP=$MOLCAS_TRAP'
  write(LuInput,'(A)') '> export MOLCAS_TRAP=ON'
  write(LuInput,'(A)') ' &RASSI &End'
  write(LuInput,'(A)') 'StOverlaps'
  write(LuInput,'(A)') 'NrOfJobIphs'
  write(LuInput,'(A)') '  2 all'
  write(LuInput,'(A)') 'IphNames'
  write(LuInput,'(A)') '  '//trim(JOB1)
  write(LuInput,'(A)') '  '//trim(JOB2)
  write(LuInput,'(A)') ' End of Input'
  if (JOB1 == 'JOBMIX') then
    write(LuInput,'(A)') '> copy $Project.JobMix JOBAUTO'
  else
    write(LuInput,'(A)') '> copy $Project.JobIph JOBAUTO'
  end if
  write(LuInput,'(A)') '> export MOLCAS_TRAP=$SL_OLD_TRAP'

else if (Request_Alaska) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Request computation of gradients.

  if (nPrint(1) >= 6) then
    call Get_cArray('Relax Method',Method,8)
    write(u6,*)
    write(u6,*) ' Slapaf requests the computation of gradients first!'
    if (iState(2) == 0) then
      write(u6,*) 'Root: ',iState(1)
      if (Method == 'MSPDFT  ') then
        NACstatesOpt(1) = iState(1)
        NACstatesOpt(2) = iState(2)
        call put_iArray('NACstatesOpt    ',NACstatesOpt,2)
        call put_lscalar('CalcNAC_Opt     ',.false.)
      end if
    else
      write(u6,*) 'Roots: ',iState(1),',',iState(2)
      if (Method == 'MSPDFT  ') then
        NACstatesOpt(1) = iState(1)
        NACstatesOpt(2) = iState(2)
        call put_iArray('NACstatesOpt    ',NACstatesOpt,2)
        ! Identify if MECI command in MC-PDFT should be used
        call put_lscalar('CalcNAC_Opt     ',.true.)
      end if
    end if
    write(u6,*)
  end if

  write(LuInput,'(A)') '>ECHO OFF'
  write(LuInput,'(A)') '>export SL_OLD_TRAP=$MOLCAS_TRAP'
  write(LuInput,'(A)') '> export MOLCAS_TRAP=ON'
  write(LuInput,'(A)') ' &Alaska &End'
  write(LuInput,'(A)') 'AUTO'
  if (iState(2) /= 0) then
    write(LuInput,'(A)') 'NAC'
    write(LuInput,'(I5,1X,I5)') iState(1),iState(2)
    write(LuInput,'(A)') 'NoCSF'
  end if
  write(LuInput,'(A)') ' End of Input'
  write(LuInput,'(A)') '> export MOLCAS_TRAP=$SL_OLD_TRAP'
  ! Repeat ECHO OFF, because it may have been turned on by alaska
  write(LuInput,'(A)') '>ECHO OFF'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Common code.

call f_inquire(Filename,Exists)
if (Exists) then
  LuSpool = 77
  LuSpool = IsFreeUnit(LuSpool)
  call Molcas_Open(LuSpool,Filename)

  do
    read(LuSpool,'(A)',iostat=istatus) Line
    if (istatus < 0) exit
    write(LuInput,'(A)') Line
  end do

  close(LuSpool)

else

  write(LuInput,'(A)') ' &Slapaf &End'
  write(LuInput,'(A)') ' End of Input'

end if
write(LuInput,'(A)') '>ECHO ON'

close(LuInput)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Start_Alaska
