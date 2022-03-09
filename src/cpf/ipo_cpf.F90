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

use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: IPOA(*)
integer(kind=iwp), intent(in) :: NVIR(*), MUL(8,8), NSYM, KLS, IFT
integer(kind=iwp) :: M, N, NSUM

NSUM = 0
do N=1,NSYM
  IPOA(N) = NSUM
  M = MUL(N,KLS)
  if (IFT < 0) then
    NSUM = NSUM+NVIR(N)*NVIR(M)
  else if (N == M) then
    NSUM = NSUM+NVIR(N)*(NVIR(N)+1)/2
  else if (N > M) then
    NSUM = NSUM+NVIR(N)*NVIR(M)
  end if
end do
IPOA(NSYM+1) = NSUM

return

end subroutine IPO_CPF
