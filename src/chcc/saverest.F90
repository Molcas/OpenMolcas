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

subroutine SaveRest(wrk,wrksize,LunAux,niter,E1old,E2old)
! this file saves 1) T1o,OE
!                 2) E1old,E2old,niter
! into RstFil file

use chcc_global, only: no, nv, PosT1o
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, LunAux, niter
real(kind=wp), intent(in) :: wrk(wrksize), E1old, E2old
integer(kind=iwp) :: len_

!open(unit=LunAux,File='RstFil',form='unformatted')
call MOLCAS_BinaryOpen_Vanilla(LunAux,'RstFil')
len_ = no*nv
call wri_chcc(LunAux,len_,wrk(PosT1o))
write(LunAux) E1old,E2old,niter
close(LunAux)

return

end subroutine SaveRest
