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

subroutine Sync_TH(TwoHam,nDens)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDens
real(kind=wp), intent(inout) :: TwoHam(nDens)
#ifdef _MOLCAS_MPP_
real(kind=wp) :: TCPU, TWall

if (.not. Is_Real_Par()) return
if (nProcs == 1) return
call BCTwoHam(TwoHam,nDens,TCPU,TWall)
#else
#include "macros.fh"
unused_var(TwoHam)
#endif

end subroutine Sync_TH
