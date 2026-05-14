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

function qRINT(N,A,C,EXPA)

use welcom, only: binom, kmax
use Constants, only: Zero, One, Half, Pi
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: qRINT
integer(kind=iwp), intent(in) :: N
real(kind=wp) :: A, C, EXPA
integer(kind=iwp) :: I, J, K, NN, nT
real(kind=wp) :: ALF, ARG, BP, Dac, F(kmax+1), FACT, FACT2, GAL, GINT, HINT, PRSUM, START, TAL

qRINT = Zero
NN = N/2+1
!write(u6,*) ' N,NN=',n,nn
BP = Half*C
START = sqrt(Pi)
PRSUM = exp(BP*BP*A+EXPA)
ALF = sqrt(A)
ARG = (BP*ALF)**2
nT = 1
F(:) = Zero
call Auxil([Arg],nT,F,nn-1)
!call RecPrt(' In qRint:Fm',' ',F,nt,nn)
GINT = Zero
Dac = -One
do I=0,N
  TAL = (-BP)**(N-I)*Binom(n,i)
  J = I/2

  if (J*2 == I) then
    Dac = Dac*Half*(real(2*J-1,kind=wp))
    FACT = ALF**(-I-1)*START*DAC
    FACT2 = BP**(I+1)*F(J+1)
    GINT = (FACT-FACT2)*TAL+GINT
  else
    GAL = One
    HINT = Zero
    do K=I-1,0,-2
      HINT = HINT+Half/A*BP**K*exp(-ARG)*GAL
      GAL = Half*real(K,kind=wp)/A*GAL
    end do
    GINT = GINT+TAL*HINT
  end if

  qRINT = GINT*PRSUM
end do

return

end function qRINT
