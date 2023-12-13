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

subroutine RandomGauss(ValMean,Sigma,iseed,nflag,buffer,Val)

use Constants, only: One, Two, Pi
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: ValMean, Sigma
real(kind=wp), intent(inout) :: buffer
real(kind=wp), intent(out) :: Val
integer(kind=iwp), intent(inout) :: iseed, nflag
real(kind=wp) :: alpha, beta, G2rad, X2pi, Z1, Z2
real(kind=wp), external :: Random_Molcas

! ValMean is the mean, and sigma is the standard deviation.
! nFlag is a binary (0,1) variable for returning the appropiate random value

! When x and y are two variables from [0, 1), uniformly
! distributed, then
!
!    cos(2*pi*x)*sqrt(-2*log(1-y))
!    sin(2*pi*x)*sqrt(-2*log(1-y))
!
! are two *independent* variables with normal distribution
! (mu = 0, sigma = 1).
! (Lambert Meertens)
! (corrected version; bug discovered by Mike Miller, fixed by LM)
if (nFlag == 0) then

  alpha = abs(Random_Molcas(iseed))
  beta = abs(Random_Molcas(iseed))

  X2pi = alpha*(Two*Pi)
  G2rad = sqrt(-Two*log(One-beta))

  Z1 = cos(X2pi)*G2rad
  Z2 = sin(X2pi)*G2rad
  Val = ValMean+Z1*Sigma
  buffer = ValMean+Z2*Sigma
  nFlag = 1

else

  Val = buffer
  nFlag = 0

end if

!write(u6,*) 'PI:' ,PPi
!write(u6,*) 'X2pi:',X2pi
!write(u6,*) 'G2grad: ',G2rad
!write(u6,*) 'buffer:',buffer
!write(u6,*) 'Z1 and Z2',Z1,Z2
!write(u6,*) 'Alpha & Beta',alpha,beta
!write(u6,*) 'VAL into RANDOM:',val

return

end subroutine RandomGauss
