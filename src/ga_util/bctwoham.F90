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
! Copyright (C) 1995, Martin Schuetz                                   *
!               1998, Roland Lindh                                     *
!***********************************************************************
!***********************************************************************
! This Module contains subroutines and functions which interface calls *
! to the Global Array Tools (GA)                                       *
!  DISTRIBUTED DATA PARALLEL VERSION for SCF                           *
!                                                                      *
! SubRoutine BCTwoHam(TwoHam,nDens,TCPU,TWall)                         *
! -> partial sum of Fock matrices                                      *
!    TwoHam: TwoHam(nDens): actual Fock matrix                         *
!----------------------------------------------------------------------*
!     written by:                                                      *
!     M. Schuetz                                                       *
!     University of Lund, Sweden, 1995                                 *
!                                                                      *
!     modified by:                                                     *
!     R. Lindh                                                         *
!     University of Lund, Sweden, 1998                                 *
!***********************************************************************

subroutine BCTwoHam(TwoHam,nDens,TCPU,TWall)

use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif

implicit none
integer(kind=iwp), intent(in) :: nDens
real(kind=wp), intent(inout) :: TwoHam(nDens)
real(kind=wp), intent(inout) :: TCPU, TWall
#ifdef _MOLCAS_MPP_
real(kind=wp) TCPU1, TWall1
real(kind=wp) TCPU2, TWall2
#include "global.fh"

if (.not. Is_Real_Par()) return
call CWTime(TCpu1,TWall1)
call GADGOP(TwoHam,nDens,'+')
call CWTime(TCpu2,TWall2)
TCPU = TCpu2-TCpu1
TWall = TWall2-TWall1

#else
#include "macros.fh"
unused_var(TwoHam)
unused_var(TCPU)
unused_var(TWall)
#endif

end subroutine BCTwoHam
