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
! Copyright (C) 2015, Luis Manuel Frutos                               *
!               2015, Ignacio Fdez. Galvan                             *
!               2015, Alessio Valentini                                *
!               2022, Isabella C. D. Merritt                           *
!               2022, Morgane Vacher                                   *
!***********************************************************************

! Modifications to calculate Wave-function Overlap using &RASSI carried out by IM
! and MV 2022

subroutine surfacehop(rc)

use Surfacehop_globals, only: firststep, Run_rassi
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: I, IAD, IDISK, istatus, ITOC15(15), LuInput, LUIPH, LuSpool, NCI, NSTATE
logical(kind=iwp) :: Exists
character(len=180) :: Line
character(len=128) :: FileName
character(len=16) :: StdIn
character(len=8) :: JOB1, JOB2, Method
real(kind=wp), allocatable :: CIBigArray(:,:)
integer(kind=iwp), external :: IsFreeUnit

#include "warnings.h"

call initial_surfacehop()
call rdinp_surfacehop()

call Get_cArray('Relax Method',Method,8)
if ((Method == 'CASPT2') .or. (Method == 'RASPT2')) then
  JOB1 = 'JOBMIX'
else
  JOB1 = 'JOBIPH'
end if

LUIPH = 20
call DANAME(LUIPH,JOB1)
IAD = 0
call IDAFILE(LUIPH,2,ITOC15,15,IAD)
call getIphInfo(LUIPH,NCI,NSTATE,ITOC15)
call mma_allocate(CIBigArray,NCI,NSTATE,Label='CIBigArray')
CIBigArray(:,:) = Zero

IDISK = ITOC15(4)
do I=1,NSTATE
  call DDAFILE(LUIPH,2,CIBigArray(:,I),NCI,IDISK)
end do
call DACLOS(LUIPH)

!call recprt('CI coefficients','',CIBigArray,NCI,NSTATE)

call tully(CIBigArray,NSTATE,NCI)

! If using CI vector product, only run tully once and quit
! Otherwise using RASSI for WF overlap
! Run_rassi is .false. if RASSI has already been run

if (.not. Run_rassi) then
  call mma_deallocate(CIBigArray)
  rc = _RC_ALL_IS_WELL_
  return
end if

! Else if tully has set Run_rassi as True, continue to run RASSI

JOB2 = 'JOBAUTO'

! If first step, cannot do overlap - just save JobIph/JobMix as JOBAUTO and return

call StdIn_Name(StdIn)
LuInput = IsFreeUnit(11)
call Molcas_Open(LuInput,StdIn)

if (firststep) then
  write(u6,*) 'First Step'
  if (JOB1 == 'JOBMIX') then
    write(u6,*) 'Saving old JobMix'
    write(LuInput,'(A)') '> copy $Project.JobMix '//trim(JOB2)
  else
    write(u6,*) 'Saving old JobIph'
    write(LuInput,'(A)') '> copy $Project.JobIph '//trim(JOB2)
  end if
  close(LuInput)
  call mma_deallocate(CIBigArray)
  rc = _RC_INVOKED_OTHER_MODULE_
  return
end if

! Otherwise, Call RASSI then rerun SURFACEHOP with same input options (check
! within tully.f90 if RASSI run yet or not)
! Call RASSI between JOB1 and JOB2

!write(u6,*) 'Calling RASSI then re-entering SURFACEHOP'

write(LuInput,'(A)') '>ECHO OFF'
write(LuInput,'(A)') '> export SH_OLD_TRAP=$MOLCAS_TRAP'
write(LuInput,'(A)') '> export MOLCAS_TRAP=ON'
write(LuInput,'(A)') '&RASSI'
write(LuInput,'(A)') 'StOverlaps'
write(LuInput,'(A)') 'NrOfJobIphs'
write(LuInput,'(A)') '  2 all'
write(LuInput,'(A)') 'IphNames'
! Note the order order of the files is opposite to e.g. SLAPAF
write(LuInput,'(A)') '  '//trim(JOB2)
write(LuInput,'(A)') '  '//trim(JOB1)
write(LuInput,'(A)') 'End of Input'
if (JOB1 == 'JOBMIX') then
  write(LuInput,'(A)') '> copy $Project.JobMix '//trim(JOB2)
else
  write(LuInput,'(A)') '> copy $Project.JobIph '//trim(JOB2)
end if
write(LuInput,'(A)') '> export MOLCAS_TRAP=$DYN_OLD_TRAP'
write(LuInput,'(A)') '>ECHO ON'

FileName = 'SURFAINP'
call f_inquire(FileName,Exists)

if (Exists) then
  LuSpool = IsFreeUnit(77)
  call Molcas_Open(LuSpool,FileName)

  do
    read(LuSpool,'(A)',iostat=istatus) Line
    if (istatus > 0) call Abend()
    if (istatus < 0) exit
    write(LuInput,'(A)') Line
  end do
  close(LuSpool)
else
  rc = _RC_INTERNAL_ERROR_
  call mma_deallocate(CIBigArray)
  close(LuInput)
  return
end if

write(LuInput,'(A)')
close(LuInput)

call mma_deallocate(CIBigArray)

rc = _RC_INVOKED_OTHER_MODULE_

end subroutine surfacehop
