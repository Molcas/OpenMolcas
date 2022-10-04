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

implicit real*8(A-H,O-Z)
real*8 FockOcc(nFockOcc)
#include "WrkSpc.fh"
#include "SysDef.fh"
character(LEN=24) Label
#ifdef _DEBUGPRINT_
#include "run_common.fh"
#endif
logical Found

! Read the generalized Fock matrix
!                                                                      *
!***********************************************************************
!                                                                      *
Label = 'FockOcc'
call qpg_dArray(Label,Found,mFockOcc)
if ((.not. Found) .or. (mFockOcc == 0)) call SysAbendMsg('get_fock_occ','Did not find:',Label)
if (mFockOcc /= nFockOcc) then
  write(6,*) 'nFockOcc=',nFockOcc
  write(6,*) 'mFockOcc=',mFockOcc
  call SysAbendMsg('get_fock_occ','mFockOcc/=nFockOcc:',Label)
end if
call Get_dArray(Label,FockOcc,nFockOcc)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
if (is_nSym == 0) then
  call get_iScalar('nSym',nSym)
  is_nSym = 1
end if
if (is_nBas == 0) then
  call Get_iArray('nBas',nBas,nSym)
  is_nBas = 1
end if
write(6,*) 'Fock occ'
ii = 1
do iIrrep=0,nSym-1
  if (nBas(iIrrep) > 0) then
    write(6,*) 'symmetry block',iIrrep
    call TriPrt(' ',' ',FockOcc(ii),nBas(iIrrep))
    ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end if
end do
#endif

return

end subroutine Get_Fock_Occ
