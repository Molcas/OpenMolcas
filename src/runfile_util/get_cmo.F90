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
! Copyright (C) Roland Lindh                                           *
!***********************************************************************
!  Get_Cmo
!
!> @brief
!>   Return the symmetry blocked MO coefficients as read from the runfile
!> @author R. Lindh
!>
!> @details
!> The utility will read the symmetry blocked MO coefficients from the runfile.
!>
!> @param[out] CMO array of symmetry blocked MO coefficients
!> @param[out] nCMO  Number of elements in the array of symmetry blocked MO coefficients
!***********************************************************************

subroutine Get_Cmo(CMO,nCMO)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nCMO
real(kind=wp) :: CMO(nCMO)
integer(kind=iwp) :: ii, iIrrep, is_nBas = 0, is_nSym = 0, mCMO, nBas(0:7) = 0, nSym = 0
character(len=24) :: Label
logical(kind=iwp) :: Found
#ifdef _DEBUGPRINT_
logical(kind=iwp), parameter :: IfTest = .true.
#else
logical(kind=iwp), parameter :: IfTest = .false.
#endif

Label = 'Last orbitals'
call qpg_dArray(Label,Found,mCmo)
if (.not. Found) call SysAbendMsg('get_CMO','Could not find',Label)
if (mCMO /= nCMO) then
  write(u6,*) 'Get_CMO: mCMO/=nCMO'
  write(u6,*) 'nCMO=',nCMO
  write(u6,*) 'mCMO=',mCMO
  call Abend()
end if
call Get_dArray(Label,CMO,nCMO)
!                                                                      *
!***********************************************************************
!                                                                      *
if (IfTest) then
  if (is_nSym == 0) then
    call Get_iScalar('nSym',nSym)
    is_nSym = 1
  end if
  if (is_nBas == 0) then
    call Get_iArray('nBas',nBas,nSym)
    is_nBas = 1
  end if
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
end if

return

end subroutine Get_Cmo
