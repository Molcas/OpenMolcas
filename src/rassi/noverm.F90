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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************

function NOVERM(N,M)

use Definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp) NOVERM
integer(kind=iwp), intent(in) :: N, M
integer(kind=iwp), save :: NOMTAB(225)
integer(kind=iwp), save :: INIT = 0
real(kind=wp) X
integer(kind=iwp) MM, IPOS, I, J, K

NOVERM = 0
if (N < 0) return
MM = M
if (2*MM > N) MM = N-M
if (MM < 0) return
NOVERM = 1
if (MM == 0) return
NOVERM = N
if (MM == 1) return
if (INIT == 0) then
  IPOS = 0
  do I=4,32
    X = real(I,kind=wp)
    do J=2,I/2
      IPOS = IPOS+1
      X = X*real(I+1-J,kind=wp)/real(J,kind=wp)
      NOMTAB(IPOS) = nint(X)
    end do
  end do
  INIT = 1
end if
if (N <= 32) then
  NOVERM = NOMTAB(((N-3)**2)/4+MM-1)
else
  X = real(NOVERM,kind=wp)
  do K=2,MM
    X = X*real(N+1-K,kind=wp)/real(K,kind=wp)
  end do
  NOVERM = nint(X)
  if (X /= real(NOVERM,kind=wp)) then
    write(u6,*) ' NOVERM: Unable to compute N over M'
    write(u6,*) ' N=',N
    write(u6,*) ' M=',M
    call ABEND()
  end if
end if

end function NOVERM
