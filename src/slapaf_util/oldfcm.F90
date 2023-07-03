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

subroutine OLDFCM(Hess,nQQ,RunOld)
!***********************************************************************
!                                                                      *
!     Object : To read in a force constant matrix from another         *
!              interphase.                                             *
!                                                                      *
!***********************************************************************

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
character*8 Method
character*(*) RunOld
logical Found
real*8, allocatable :: Hess(:)

! Prologue

! Set runfile to be the one according to the character string RUNOLD
call NameRun(RunOld)

! Get the method used in the old calculation
call Get_cArray('Relax Method',Method,8)

! Get the final energy obtained by the last calculation
!call Get_Energy(Energy)
call Get_dScalar('Last energy',Energy)

! Get the number of internal coordinates
call Get_iScalar('No of Internal coordinates',iInter)
if (iInter <= 0) then
  call WarningMessage(2,'OldFCM: iInter <= 0')
  write(6,*) 'iInter=',iInter
  call Abend()
end if

! Get the force constant matrix
call qpg_dArray('Hess',Found,nHess)
if ((.not. Found) .or. (nHess == 0)) then
  call SysAbendmsg('OldFcm','Did not find:','Hess')
end if

call mma_Allocate(Hess,nHess,Label='Hess')
call get_dArray('Hess',Hess,nHess)

lHess = iInter**2
if (nHess /= lHess) then
  call WarningMessage(2,'OldFCM: nHess /= lHess')
  write(6,*) 'nHess,lHess=',nHess,lHess
  call Abend()
end if

! Reset runfile to be RUNFILE
call NameRun('#Pop')

! Echo the input information
#ifdef _DEBUGPRINT_
write(6,*)
write(6,'(6X,A)') 'SLAPAF has been supplied with an old force constant matrix.'
write(6,'(6X,3A)') 'It is based on ',Method,' calculations.'
write(6,'(6X,A,F18.10)') 'The final energy was',Energy
call RecPrt(' OldFcm',' ',Hess,iInter,iInter)
#endif

nQQ = iINter

! Epilogue, end
return

end subroutine OLDFCM
