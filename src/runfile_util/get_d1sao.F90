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

subroutine Get_D1sao(D1sao,nD1sao)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nD1sao
real(kind=wp) :: D1sao(nD1sao)
integer(kind=iwp) :: nDens
logical(kind=iwp) :: Found
character(len=24) :: Label
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nBas(0:7) = -1, nSym = -1
#endif

! Read the variational 1st order density matrix
! density matrix in AO/SO basis
!                                                                      *
!***********************************************************************
!                                                                      *
Label = 'D1sao'
call qpg_dArray(Label,Found,nDens)
if ((.not. Found) .or. (nDens == 0)) call SysAbendMsg('get_d1sao','Did not find',Label)
if (nDens /= nD1sao) then
  write(u6,*) 'Get_D1sao: nDens/=nD1sao'
  call Abend()
end if
call Get_dArray(Label,D1sao,nD1sao)
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
    call TriPrt(' ',' ',D1sao(ii),nBas(iIrrep))
    ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end if
end do
#endif

return

end subroutine Get_D1sao
