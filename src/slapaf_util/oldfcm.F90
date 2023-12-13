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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine OLDFCM(Hess,nQQ,RunOld)
!***********************************************************************
!                                                                      *
!     Object : To read in a force constant matrix from another         *
!              interphase.                                             *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), allocatable, intent(out) :: Hess(:)
integer(kind=iwp), intent(out) :: nQQ
character(len=*), intent(in) :: RunOld
integer(kind=iwp) :: iInter, lHess, nHess
real(kind=wp) :: Energy
character(len=8) :: Method
logical(kind=iwp) :: Found

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
  write(u6,*) 'iInter=',iInter
  call Abend()
end if

! Get the force constant matrix
call qpg_dArray('Hess',Found,nHess)
if ((.not. Found) .or. (nHess == 0)) call SysAbendmsg('OldFcm','Did not find:','Hess')

call mma_Allocate(Hess,nHess,Label='Hess')
call get_dArray('Hess',Hess,nHess)

lHess = iInter**2
if (nHess /= lHess) then
  call WarningMessage(2,'OldFCM: nHess /= lHess')
  write(u6,*) 'nHess,lHess=',nHess,lHess
  call Abend()
end if

! Reset runfile to be RUNFILE
call NameRun('#Pop')

! Echo the input information
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,'(6X,A)') 'SLAPAF has been supplied with an old force constant matrix.'
write(u6,'(6X,3A)') 'It is based on ',Method,' calculations.'
write(u6,'(6X,A,F18.10)') 'The final energy was',Energy
call RecPrt(' OldFcm',' ',Hess,iInter,iInter)
#endif

nQQ = iINter

! Epilogue, end
return

end subroutine OLDFCM
