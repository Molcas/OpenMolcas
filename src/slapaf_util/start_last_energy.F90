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

subroutine Start_Last_Energy()

use RunFile_procedures, only: Get_Coord_New
use stdalloc, only: mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: iGO, lengthlast, LuInput, nCoord, nSaddle
logical(kind=iwp) :: FoundLastEn, Saddle
character(len=16) :: StdIn
character(len=8) :: Method
real(kind=wp), allocatable :: CN(:,:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
if (nPrint(1) >= 6) then
  write(u6,*)
  write(u6,*) ' Slapaf requests the last energy to be computed!'
  write(u6,*)
end if

LuInput = 11
LuInput = IsFreeUnit(LuInput)
call StdIn_Name(StdIn)
call Molcas_Open(LuInput,StdIn)

call Get_iScalar('Grad ready',iGO)

write(LuInput,'(A)') '>ECHO OFF'
write(LuInput,'(A)') '>export SL_OLD_TRAP=$MOLCAS_TRAP'
write(LuInput,'(A)') '>export MOLCAS_TRAP=ON'

call qpg_dArray('Saddle',Saddle,nSaddle)
if (Saddle) then

  ! Clean up among the runfiles, etc.

  write(LuInput,'(A)') '>COPY $Project$SubProject.RunFile $Project.RunFile'

  call qpg_cArray('LastEnergyMethod',FoundLastEn,lengthlast)
  if (FoundLastEn) then
    call Get_cArray('LastEnergyMethod',Method,8)
  else
    call Get_cArray('Relax Method',Method,8)
  end if

  if ((Method == 'CASSCF') .or. (Method == 'RASSCF') .or. (Method == 'CASSCFSA') .or. (Method == 'RASSCFSA') .or. &
      (Method == 'CASPT2') .or. (Method == 'RASPT2')) then

    write(LuInput,'(A)') '>COPY $Project$SubProject.JobIph $Project.JobIph'
    write(LuInput,'(A)') '>RM $Project.Reac.JobIph'
    write(LuInput,'(A)') '>RM $Project.Prod.JobIph'
  end if
  write(LuInput,'(A)') '>RM $Project.Reac.RunFile'
  write(LuInput,'(A)') '>RM $Project.Prod.RunFile'
  write(LuInput,'(A)') '>export SubProject='
  write(LuInput,'(A)') '>export MOLCAS_SADDLE=0'

  ! Put the final structure as the reference structure

  call Get_Coord_New(CN,nCoord)
  call NameRun('RUNREAC')
  call Put_dArray('Ref_Geom',CN,3*nCoord)
  call NameRun('#Pop')
  call NameRun('RUNPROD')
  call Put_dArray('Ref_Geom',CN,3*nCoord)
  call NameRun('#Pop')
  call mma_deallocate(CN)
end if

if (iGo <= 1) then

  ! Normal behaviour

  write(LuInput,'(A)') ' &Last_Energy &End'
  write(LuInput,'(A)') 'End of Input'

else

  ! A CI ISC search has been performed. For the moment I care only
  ! for the ISC. Once the rest is done I will take care of that as well.

  write(LuInput,'(A)') '>COPY $OldProject.Seward.Input State1.Seward.Input'
  write(LuInput,'(A)') '>COPY $OldProject.Seward.Input State2.Seward.Input'
  write(LuInput,'(A)') '>COPY $OldProject.RunFile State1.RunFile'
  write(LuInput,'(A)') '>COPY $OldProject.RunFile State2.RunFile'
  write(LuInput,'(A)') '>RM molcas.env'
  write(LuInput,'(A)') '>export Project=State1'
  write(LuInput,'(A)') ' &Last_Energy &End'
  write(LuInput,'(A)') 'End of Input'
  write(LuInput,'(A)') '>RM molcas.env'
  write(LuInput,'(A)') '>export Project=State2'
  write(LuInput,'(A)') ' &Last_Energy &End'
  write(LuInput,'(A)') 'End of Input'
end if
write(LuInput,'(A)') '>export MOLCAS_TRAP=$SL_OLD_TRAP'

write(LuInput,'(A)') '>ECHO ON'

close(LuInput)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Start_Last_Energy
