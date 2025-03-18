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

integer function iBion(M,N)
! Compute binomial coefficients recursively

real*8 Temp1, Temp2, Temp3, Temp4

if ((N < 0) .or. (M < N)) then
  write(6,*) 'Wrong params is iBion',M,N
  call Abend()
end if
if ((N == 0) .or. (M == 0)) then
  iBion = 1
  return
else if (N == 1) then
  iBion = M
  return
else if (N == 2) then
  iBion = M*(M-1)/2
  return
end if
K1 = max((M-N),N)
K2 = min((M-N),N)
Temp1 = 1.0d0
do K=1,K2
  Temp2 = dble(K1+K)
  Temp3 = dble(K)
  Temp4 = Temp2/Temp3
  Temp1 = Temp1*Temp4
end do
iBion = nint(Temp1)

return

end function iBion
