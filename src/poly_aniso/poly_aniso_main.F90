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

subroutine poly_aniso_open(iReturn)

implicit none
#include "warnings.h"
integer :: nneq, exch, neqv, nmax, nLoc, nCenter, nT, nH, nTempMagn, nDir, nDirZee, nMult, nPair, MxRank1, MxRank2
integer :: iReturn
logical :: dbg, old_aniso_format

iReturn = 0
dbg = .false.

write(6,'(A)') 'POLY_ANISO (OPEN):'
write(6,'(A)') 'by:   Liviu Unugur       (chmlu@nus.edu.sg)'
write(6,'(A)') 'and   Liviu F. Chibotaru (Liviu.Chibotaru@kuleuven.be)'
write(6,'(A)') 'Last updated - 2 July 2018'

! find if we are using old format or the new one
if (dbg) write(6,*) 'Enter find_aniso_format'
call find_aniso_format(old_aniso_format)
if (dbg) write(6,*) 'Exit find_aniso_format'

! initialize some important variables
if (dbg) write(6,*) 'Enter fetch_init_const'
call fetch_init_const(nneq,neqv,nmax,exch,nLoc,nCenter,nT,nH,nTempMagn,nDir,nDirZee,nMult,nPair,MxRank1,MxRank2,old_aniso_format, &
                      iReturn)
if (dbg) write(6,*) 'Exit fetch_init_const'

if (dbg) write(6,*) 'pa_open::             nneq=',nneq
if (dbg) write(6,*) 'pa_open::             neqv=',neqv
if (dbg) write(6,*) 'pa_open::             nmax=',nmax
if (dbg) write(6,*) 'pa_open::             exch=',exch
if (dbg) write(6,*) 'pa_open::             nLoc=',nLoc
if (dbg) write(6,*) 'pa_open::          nCenter=',nCenter
if (dbg) write(6,*) 'pa_open::               nT=',nT
if (dbg) write(6,*) 'pa_open::               nH=',nH
if (dbg) write(6,*) 'pa_open::        nTempMagn=',nTempMagn
if (dbg) write(6,*) 'pa_open::             nDir=',nDir
if (dbg) write(6,*) 'pa_open::          nDirZee=',nDirZee
if (dbg) write(6,*) 'pa_open::            nMult=',nMult
if (dbg) write(6,*) 'pa_open::            nPair=',nPair
if (dbg) write(6,*) 'pa_open::          MxRank1=',MxRank1
if (dbg) write(6,*) 'pa_open::          MxRank2=',MxRank2
if (dbg) write(6,*) 'pa_open::          iReturn=',iReturn
if (dbg) write(6,*) 'pa_open:: old_aniso_format=',old_aniso_format
if (dbg) call xFlush(6)

if (iReturn /= 0) then
  write(6,*) 'ERROR: something went wrong during initialization of main variables'
  write(6,*) 'Have to quit now'
  call Quit(_RC_GENERAL_ERROR_)
end if

! main program
if (dbg) write(6,*) 'Enter poly_aniso_1'
if (dbg) call xFlush(6)
call POLY_ANISO_1(nneq,neqv,nmax,exch,nLoc,nCenter,nT,nH,nTempMagn,nDir,nDirZee,nMult,nPair,MxRank1,MxRank2,old_aniso_format, &
                  iReturn)
if (dbg) write(6,*) 'Exit poly_aniso_1'

if (iReturn /= 0) then
  write(6,*) 'ERROR: something went wrong during the execution of the POLY_ANISO_1'
  call Quit(_RC_GENERAL_ERROR_)
end if
if (dbg) call xFlush(6)

return

end subroutine poly_aniso_open
