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

subroutine GetRest(wrk,wrksize,LunAux,niter,E1old,E2old)
! this file reads 1) T1o
!                 2) E1old,E2old,niter
! from RstFil file

use chcc_global, only: no, nv, PosT1o
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, LunAux
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
integer(kind=iwp), intent(out) :: niter
real(kind=wp), intent(out) :: E1old, E2old
integer(kind=iwp) :: len_

!open(unit=LunAux,File='RstFil',form='unformatted')
call MOLCAS_BinaryOpen_Vanilla(LunAux,'RstFil')
len_ = no*nv
call rea_chcc(LunAux,len_,wrk(PosT1o))
read(LunAux) E1old,E2old,niter
close(LunAux)

return

end subroutine GetRest
