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
!               2022, Isabella C. D. Merritt                            *
!               2022, Morgane Vacher                                   *
!***********************************************************************

! Modifications to calculate Wave-function Overlap using &RASSI carried out by IM
! and MV 2022

subroutine surfacehop(rc)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6
use Tully_variables, only: rassi_ovlp, firststep, Run_rassi

implicit none
integer(kind=iwp), intent(out) :: rc
#include "warnings.h"

integer(kind=iwp) :: NSTATE, LUIPH, IAD, ITOC15(15), NCI, IDISK, I
real(kind=wp), allocatable :: CIBigArray(:)

character(len=16) :: StdIn
character(len=128) :: FileName
character(len=180) :: Line
integer(kind=iwp) :: LuInput, LuSpool, istatus
integer(kind=iwp), external :: IsFreeUnit
logical(kind=iwp) :: Exists

call initial_surfacehop()
call rdinp_surfacehop()

LUIPH = 20
call DANAME(LUIPH,'JOBIPH')
IAD = 0
call IDAFILE(LUIPH,2,ITOC15,15,IAD)
call getIphInfo(LUIPH,NCI,NSTATE,ITOC15)
call mma_allocate(CIBigArray,NCI*NSTATE)
CIBigArray(:) = Zero

IDISK = ITOC15(4)
do I=1,NSTATE
  call DDAFILE(LUIPH,2,CIBigArray(NCI*(I-1)+1),NCI,IDISK)
end do
call DACLOS(LUIPH)

!call recprt('CI coefficients','',CIBigArray,NCI,NSTATE)


! If using CI vector product, only run tully once and quit
if (.not. rassi_ovlp) then
  call tully(CIBigArray,NSTATE,NCI)
  call mma_deallocate(CIBigArray)
  rc = _RC_ALL_IS_WELL_
  return
end if

! Otherwise using RASSI for WF overlap

call tully(CIBigArray,NSTATE,NCI)

if (.not. Run_rassi) then ! RASSI Already Run
  call mma_deallocate(CIBigArray)
  rc = _RC_ALL_IS_WELL_
  return
end if

! Else if tully has set Run_rassi as True, continue to run RASSI

! If first step, cannot do overlap - just save JobIph as JobOld and return

if (firststep) then
  write(u6,*) 'First Step'
  LuInput=11
  LuInput=IsFreeUnit(LuInput)
  write(u6,*) 'Saving old JobIPH'
  Call StdIn_Name(StdIn)
  Call Molcas_Open(LuInput,StdIn)
  Write (LuInput,'(A)') ' >copy $Project.JobIph $Project.JobIph.Old'
  Close(LuInput)
  call mma_deallocate(CIBigArray)
  rc=_RC_INVOKED_OTHER_MODULE_
  return
end if

! Otherwise, Call RASSI then rerun SURFACEHOP with same input options (check
! within tully.f90 if RASSI run yet or not)
! Call RASSI between .JobIph and .JobOld

LuInput=11
LuInput=IsFreeUnit(LuInput)
!write(u6,*) 'Calling RASSI then re-entering SURFACEHOP'

Call StdIn_Name(StdIn)
Call Molcas_Open(LuInput,StdIn)

Write (LuInput,'(A)') '>copy $Project.JobIph.Old JOB001'
Write (LuInput,'(A)') '>copy $Project.JobIph JOB002'

Write (LuInput,'(A)') '&RASSI &End'
Write (LuInput,'(A)') 'Nr of JobIPhs'
Write (LuInput,'(A)') '2 all'
Write (LuInput,'(A)') 'STOV'
Write (LuInput,'(A)') 'End of Input'
Write (LuInput,'(A)') '> copy $Project.JobIph $Project.JobIph.Old'

FileName = 'SURFAINP'
call f_inquire(FileName,Exists)

if (Exists) then
  LuSpool = 77
  LuSpool = IsFreeUnit(LuSpool)
  call Molcas_Open(LuSpool, FileName)

  do
    read(LuSpool,'(A)',iostat=istatus) Line
    if (istatus > 0) call Abend()
    if (istatus < 0) exit
    write(LuInput,'(A)') Line
  end do
  Close(LuSpool)
else
  rc = _RC_INTERNAL_ERROR_
  call mma_deallocate(CIBigArray)
  Close(LuInput)
  return
end if

Write (LuInput,'(A)') ''
Close(LuInput)

call mma_deallocate(CIBigArray)

call Finish(_RC_INVOKED_OTHER_MODULE_)

return

end subroutine surfacehop
