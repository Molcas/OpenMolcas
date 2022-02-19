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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine IPO_CPF(IPOA,NVIR,MUL,NSYM,KLS,IFT)

implicit real*8(A-H,O-Z)
dimension IPOA(*), NVIR(*), MUL(8,8)

NSUM = 0
do N=1,NSYM
  IPOA(N) = NSUM
  M = MUL(N,KLS)
  if (IFT >= 0) GO TO 20
  NSUM = NSUM+NVIR(N)*NVIR(M)
  GO TO 10
20 if (N-M < 0) then
    GO TO 10
  else if (N-M == 0) then
    GO TO 11
  else
    GO TO 12
  end if
11 NSUM = NSUM+NVIR(N)*(NVIR(N)+1)/2
  GO TO 10
12 NSUM = NSUM+NVIR(N)*NVIR(M)
10 continue
end do
IPOA(NSYM+1) = NSUM

return

end subroutine IPO_CPF
