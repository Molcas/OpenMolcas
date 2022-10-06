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

subroutine Get_Cmo_(CMO,nCMO)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nCMO
real(kind=wp) :: CMO(nCMO)
integer(kind=iwp) :: mCMO
character(len=24) :: Label
logical(kind=iwp) :: Found
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nBas(0:7) = -1, nSym = -1
#endif

Label = 'Last orbitals'
call qpg_dArray(Label,Found,mCmo)
if (.not. Found) then
  Label = 'Guessorb'
  call qpg_dArray(Label,Found,mCmo)
  if (.not. Found) call SysAbendMsg('get_CMO','Could not find',Label)
end if
if (mCMO /= nCMO) then
  write(u6,*) 'Get_CMO_: mCMO/=nCMO'
  write(u6,*) 'nCMO=',nCMO
  write(u6,*) 'mCMO=',mCMO
  call Abend()
end if

call Get_dArray(Label,CMO,nCMO)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
if (nSym < 0) call Get_iScalar('nSym',nSym)
if (nBas(0) < 0) call Get_iArray('nBas',nBas,nSym)
write(u6,*) ' Input Orbitals from RUNFILE'
write(u6,*)
ii = 1
do iIrrep=0,nSym-1
  if (nBas(iIrrep) > 0) then
    write(u6,*) ' Symmetry Block',iIrrep
    call RecPrt(' ',' ',CMO(ii),nBas(iIrrep),nBas(iIrrep))
    write(u6,*)
  end if
  ii = ii+nBas(iIrrep)**2
end do
#endif

return

end subroutine Get_Cmo_
