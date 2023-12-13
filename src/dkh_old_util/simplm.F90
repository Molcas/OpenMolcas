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

function SIMPLM(NMP,F,X)
! Calculates  Intg/xmin,xmax/(F(x) dx)
! if the F function has been evaluated in a logarithmic mesh.
!
! Note that the integral equals Ingt/xmin,xmax/(F.x d(ln x)),
! which is what is actually calculated.
!
! NMP:    number of mesh points
! F:      array with the function values
! X:      array with the variable values of the log. mesh

use Constants, only: Zero, Three, Four, Six, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: SIMPLM
integer(kind=iwp), intent(in) :: NMP
real(kind=wp), intent(in) :: F(NMP), X(NMP)
integer(kind=iwp) :: I, IP1, IP2, J, N, NM2
real(kind=wp) :: DELTA, DLM, SUM_

! checks the precision of the incr. of the log.mesh
DLM = log(X(2))-log(X(1))
do J=2,5
  DELTA = log(X(J+1))-log(X(J))
  if (abs(DELTA-DLM) >= 1.0e-8_wp) then
    write(u6,602)
    call Abend()
  end if
end do

if (mod(NMP,2) /= 0) then
  N = NMP
else
  N = NMP-1
end if
! Odd number of points
SUM_ = Zero
NM2 = N-2
do I=1,NM2,2
  IP1 = I+1
  IP2 = I+2
  SUM_ = SUM_+F(I)*X(I)+Four*F(IP1)*X(IP1)+F(IP2)*X(IP2)
end do
SIMPLM = SUM_*DLM/Three
if (N == NMP) return

! even number of points: add up the remaining
SUM_ = 2.5_wp*F(NMP)*X(NMP)+Four*F(N)*X(N)-Half*F(N-1)*X(N-2)
SIMPLM = SIMPLM+SUM_*DLM/Six

return

602 format(' SIMPLM: Increment of the log mesh not constant')

end function SIMPLM
