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

subroutine Get_Fock_Occ(FockOcc,nFockOcc)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nFockOcc
real(kind=wp) :: FockOcc(nFockOcc)
integer(kind=iwp) :: mFockOcc
logical(kind=iwp) :: Found
character(len=24) :: Label
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nBas(0:7) = -1, nSym = -1
#endif

! Read the generalized Fock matrix
!                                                                      *
!***********************************************************************
!                                                                      *
Label = 'FockOcc'
call qpg_dArray(Label,Found,mFockOcc)
if ((.not. Found) .or. (mFockOcc == 0)) call SysAbendMsg('get_fock_occ','Did not find:',Label)
if (mFockOcc /= nFockOcc) then
  write(u6,*) 'nFockOcc=',nFockOcc
  write(u6,*) 'mFockOcc=',mFockOcc
  call SysAbendMsg('get_fock_occ','mFockOcc/=nFockOcc:',Label)
end if
call Get_dArray(Label,FockOcc,nFockOcc)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
if (nSym < 0) call get_iScalar('nSym',nSym)
if (nBas(0) < 0) call Get_iArray('nBas',nBas,nSym)
write(u6,*) 'Fock occ'
ii = 1
do iIrrep=0,nSym-1
  if (nBas(iIrrep) > 0) then
    write(u6,*) 'symmetry block',iIrrep
    call TriPrt(' ',' ',FockOcc(ii),nBas(iIrrep))
    ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end if
end do
#endif

return

end subroutine Get_Fock_Occ
