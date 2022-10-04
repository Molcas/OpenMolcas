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

subroutine Get_D1ao_Var(D1ao,nD1ao)

implicit real*8(A-H,O-Z)
#include "WrkSpc.fh"
#include "SysDef.fh"
character*24 Label
#ifdef _DEBUGPRINT_
#include "run_common.fh"
#endif
logical Found
integer nD1ao
real*8 D1ao(nD1ao)

! Read the variational 1st order density matrix
! density matrix in AO/SO basis
!                                                                      *
!***********************************************************************
!                                                                      *
Label = 'D1aoVar'
call qpg_dArray(Label,Found,nDens)
if ((.not. Found) .or. (nDens == 0)) then
  call Get_D1ao(D1ao,nD1ao)
else if (nDens == nD1ao) then
  call get_dArray(Label,D1ao,nD1ao)
else
  write(6,*) 'Get_D1ao_Var: nDens/=nD1ao'
  write(6,*) 'nDens=',nDens
  write(6,*) 'nD1ao=',nD1ao
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
if (is_nSym == 0) then
  call Get_iScalar('nSym',nSym)
  is_nSym = 1
end if
if (is_nBas == 0) then
  call Get_iArray('nBas',nBas,nSym)
  is_nBas = 1
end if
write(6,*) 'variational 1st order density matrix'
ii = 1
do iIrrep=0,nSym-1
  if (nBas(iIrrep) > 0) then
    write(6,*) 'symmetry block',iIrrep
    call TriPrt(' ',' ',D1ao(ii),nBas(iIrrep))
    ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end if
end do
#endif

return

end subroutine Get_D1ao_Var
