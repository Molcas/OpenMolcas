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

subroutine Get_D1ao_ab(ipD1ao,nDens)

implicit real*8(A-H,O-Z)
#include "WrkSpc.fh"
#include "SysDef.fh"
character*24 Label
logical Found

! Read the variational 1st order density matrix
! density matrix in AO/SO basis
!                                                                      *
!***********************************************************************
!                                                                      *
Label = 'D1ao_ab'
call qpg_dArray(Label,Found,nDens)
if ((.not. Found) .or. (nDens == 0)) call SysAbendMsg('get_d1ao_ab','Could not locate:',Label)
call GetMem('Dens_ab','Allo','Real',ipD1ao,nDens)
call get_dArray(Label,Work(ipD1ao),nDens)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Get_D1ao_ab
