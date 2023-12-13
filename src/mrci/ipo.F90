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

subroutine IPO(IPOA,NVIR,MUL,NSYM,KLS,IFT)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NVIR(*), MUL(8,8), NSYM, KLS, IFT
integer(kind=iwp), intent(out) :: IPOA(NSYM+1)
integer(kind=iwp) :: M, N, NSUM

NSUM = 0
if (IFT < 0) then
  do N=1,NSYM
    IPOA(N) = NSUM
    M = MUL(N,KLS)
    NSUM = NSUM+NVIR(N)*NVIR(M)
  end do
else
  if (KLS == 1) then
    do N=1,NSYM
      IPOA(N) = NSUM
      NSUM = NSUM+(NVIR(N)*(NVIR(N)+1))/2
    end do
  else
    do N=1,NSYM
      IPOA(N) = NSUM
      M = MUL(N,KLS)
      if (N > M) then
        NSUM = NSUM+NVIR(N)*NVIR(M)
      end if
    end do
  end if
end if
IPOA(NSYM+1) = NSUM

return

end subroutine IPO
