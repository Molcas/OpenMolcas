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

integer function IBION_LUCIA(M,N)

implicit none
integer M, N
integer, external :: IBINOM

! PAM05:
! The following code does not always work.
! Replaced by call to my 'NOVERM' routine, renamed IBINOM.
!
! BIONOMIAL COEFFICIENT (M / N ) = IFAC(M)/(IFAC(M-N)*IFAC(N))
!
!IB = 1
!if (M-N >= N) then
!  do K=(M-N+1),M
!    IB = IB*K
!  end do
!  IB = IB/IFAC(N)
!else
!  do K=N+1,M
!    IB = IB*K
!  end do
!  IB = IB/IFAC(M-N)
!end if
!
!IBION_LUCIA = IB

IBION_LUCIA = IBINOM(M,N)

end function IBION_LUCIA
