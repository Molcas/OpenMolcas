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
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************

subroutine Diagr(A,N,EIG,EW,SINV,AUX,AUXI)

use Constants, only: One, Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: A(*), SINV(N,N)
real(kind=wp), intent(out) :: EIG(N,N), AUX(N,N)
real(kind=wp), intent(_OUT_) :: EW(*), AUXI(*)

if (n == 0) return
call Square(A,Aux,n,1,n)
call DGEMM_('N','N',n,n,n,One,Aux,n,Sinv,n,Zero,Eig,n)
call dGemm_tri('T','N',n,n,n,One,SINV,n,EIG,n,Zero,AUXI,n)

call unitmat(Eig,N)
call dCopy_(N*(N+1)/2,AUXI,1,AUX,1)
!call NIDiag(AUXI,EIG,N,N)
call NIDiag_New(AUXI,EIG,N,N)
call vEig(N,AUXI,EW)
call SortEig(EW,Eig,n,n,1,.false.)

return

end subroutine Diagr
