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

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nD1ao
real(kind=wp), intent(out) :: D1ao(nD1ao)
integer(kind=iwp) :: nDens
logical(kind=iwp) :: Found
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ii, iIrrep, nBas(0:7) = -1, nSym = -1
#endif
character(len=*), parameter :: Label = 'D1aoVar'

! Read the variational 1st order density matrix
! density matrix in AO/SO basis
!                                                                      *
!***********************************************************************
!                                                                      *
call qpg_dArray(Label,Found,nDens)
if ((.not. Found) .or. (nDens == 0)) then
  call Get_dArray_chk('D1ao',D1ao,nD1ao)
else if (nDens == nD1ao) then
  call get_dArray(Label,D1ao,nD1ao)
else
  write(u6,*) 'Get_D1ao_Var: nDens/=nD1ao'
  write(u6,*) 'nDens=',nDens
  write(u6,*) 'nD1ao=',nD1ao
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
if (nSym < 0) call Get_iScalar('nSym',nSym)
if (nBas(0) < 0) call Get_iArray('nBas',nBas,nSym)
write(u6,*) 'variational 1st order density matrix'
ii = 1
do iIrrep=0,nSym-1
  if (nBas(iIrrep) > 0) then
    write(u6,*) 'symmetry block',iIrrep
    call TriPrt(' ',' ',D1ao(ii),nBas(iIrrep))
    ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end if
end do
#endif

return

end subroutine Get_D1ao_Var
