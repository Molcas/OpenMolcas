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

subroutine CpLabr(A,B,L,M,N,IA,IB,C,IC,IER)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L, M, N, IA, IB, IC
real(kind=wp), intent(in) :: A(IA,M), B(IB,N)
real(kind=wp), intent(inout) :: C(IC,N)
integer(kind=iwp), intent(out) :: IER

#ifdef _DEBUGPRINT_
if ((IA >= L) .and. (IB >= M) .and. (IC >= L)) then
#endif
  IER = 0
  call DGEMM_('N','N',L,M,N,One,A,IA,B,IB,One,C,IC)
#ifdef _DEBUGPRINT_
else
  IER = 129
end if
#endif

return

end subroutine CpLabr
