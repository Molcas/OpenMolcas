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

subroutine GR_DArray(Array,nArray)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArray
real(kind=wp), intent(inout) :: Array(nArray)
#ifdef _MOLCAS_MPP_
real(kind=wp) :: TCpu1, TCpu2, TWall1, TWall2

if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
call CWTime(TCpu1,TWall1)
call GADGOP(Array,nArray,'+')
call CWTime(TCpu2,TWall2)
#else
#include "macros.fh"
unused_var(Array)
#endif

end subroutine GR_DArray
