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

subroutine ClsFls_MCLR()
!***********************************************************************
!                                                                      *
!     Open files.                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use MCLR_Data, only: FnMck, LuCSF2SD, LuJob, LuMck, LuQDat, LuTemp, LuTri1, SA
use input_mclr, only: iMethod, RASSI, TwoStep
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iOpt, iRC
logical(kind=iwp) :: DoCholesky
integer(kind=iwp), external :: AixRm

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
if (iMethod == 2) then
  call DaClos(LuCSF2sd)
  !---  close the JOBIPH file -----------------------------------------*
  call DaClos(LuJob)
end if
call DaClos(LuTemp)
!---  close the ORDINT file -------------------------------------------*
call DecideonCholesky(DoCholesky)
if (.not. DoCholesky) then
  iRc = -1
  call ClsOrd(iRc)
  if (iRc /= 0) then
    write(u6,*) 'ClsFls: Error closing ORDINT'
    call Abend()
  end if
end if
call DaClos(LuTri1)
if (TwoStep) then
  call DaClos(LuQDAT)
  !call DaClos(LuMOTRA)
end if

! Close the MckInt file or Remove the MCKINT file if SA----------------*
! Do not remove file if we are producing data on the MckInt file for
! the RASSI module!

if (SA .and. (.not. RASSI)) then
  ! What the...? No control at all on what file is being removed!
  !call DaEras(LuMck)
  call DaClos(LuMck)
  iRC = AixRM(FnMck)
else
  iRc = -1
  iOpt = 0
  call ClsMck(iRc,iOpt)
  if (iRc /= 0) then
    write(u6,*) 'ClsFls: Error closing MCKINT'
    call Abend()
  end if
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine ClsFls_MCLR
