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

use Definitions, only: iwp, u6

implicit none
#include "warnings.h"
integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: exch, MxRank1, MxRank2, nCenter, nDir, nDirZee, neqv, nH, nLoc, nmax, nMult, nneq, nPair, nT, nTempMagn
logical(kind=iwp) :: old_aniso_format

iReturn = 0

write(u6,'(A)') 'POLY_ANISO (OPEN):'
write(u6,'(A)') 'by:   Liviu Unugur       (chmlu@nus.edu.sg)'
write(u6,'(A)') 'and   Liviu F. Chibotaru (Liviu.Chibotaru@kuleuven.be)'
write(u6,'(A)') 'Last updated - 2 July 2018'

! find if we are using old format or the new one
#ifdef _DEBUGPRINT_
write(u6,*) 'Enter find_aniso_format'
#endif
call find_aniso_format(old_aniso_format)
#ifdef _DEBUGPRINT_
write(u6,*) 'Exit find_aniso_format'
#endif

! initialize some important variables
#ifdef _DEBUGPRINT_
write(u6,*) 'Enter fetch_init_const'
#endif
call fetch_init_const(nneq,neqv,nmax,exch,nLoc,nCenter,nT,nH,nTempMagn,nDir,nDirZee,nMult,nPair,MxRank1,MxRank2,old_aniso_format, &
                      iReturn)
#ifdef _DEBUGPRINT_
write(u6,*) 'Exit fetch_init_const'

write(u6,*) 'pa_open::             nneq=',nneq
write(u6,*) 'pa_open::             neqv=',neqv
write(u6,*) 'pa_open::             nmax=',nmax
write(u6,*) 'pa_open::             exch=',exch
write(u6,*) 'pa_open::             nLoc=',nLoc
write(u6,*) 'pa_open::          nCenter=',nCenter
write(u6,*) 'pa_open::               nT=',nT
write(u6,*) 'pa_open::               nH=',nH
write(u6,*) 'pa_open::        nTempMagn=',nTempMagn
write(u6,*) 'pa_open::             nDir=',nDir
write(u6,*) 'pa_open::          nDirZee=',nDirZee
write(u6,*) 'pa_open::            nMult=',nMult
write(u6,*) 'pa_open::            nPair=',nPair
write(u6,*) 'pa_open::          MxRank1=',MxRank1
write(u6,*) 'pa_open::          MxRank2=',MxRank2
write(u6,*) 'pa_open::          iReturn=',iReturn
write(u6,*) 'pa_open:: old_aniso_format=',old_aniso_format
call xFlush(u6)
#endif

if (iReturn /= 0) then
  write(u6,*) 'ERROR: something went wrong during initialization of main variables'
  write(u6,*) 'Have to quit now'
  call Quit(_RC_GENERAL_ERROR_)
end if

! main program
#ifdef _DEBUGPRINT_
write(u6,*) 'Enter poly_aniso_1'
call xFlush(u6)
#endif
call POLY_ANISO_1(nneq,neqv,nmax,exch,nLoc,nCenter,nT,nH,nTempMagn,nDir,nDirZee,nMult,nPair,MxRank1,MxRank2,old_aniso_format, &
                  iReturn)
#ifdef _DEBUGPRINT_
write(u6,*) 'Exit poly_aniso_1'
#endif

if (iReturn /= 0) then
  write(u6,*) 'ERROR: something went wrong during the execution of the POLY_ANISO_1'
  call Quit(_RC_GENERAL_ERROR_)
end if
#ifdef _DEBUGPRINT_
call xFlush(u6)
#endif

return

end subroutine poly_aniso_open
