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

implicit real*8(a-h,o-z)
real*8 A(IA,M), B(IB,N), C(IC,N)

#ifdef _DEBUGPRINT_
if ((IA >= L) .and. (IB >= M) .and. (IC >= L)) GO TO 5
IER = 129
GO TO 9000
5 continue
#endif
IER = 0
call DGEMM_('N','N',L,M,N,1.0d0,A,IA,B,IB,1.0d0,C,IC)
#ifdef _DEBUGPRINT_
9000 continue
#endif

return

end subroutine CpLabr
