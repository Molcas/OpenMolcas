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

subroutine GramSchmidt(S,C,nDim,itype,Center,iRestrict)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, itype(nDim), Center(nDim), iRestrict
real(kind=wp), intent(inout) :: S(nDim,nDim), C(nDim,nDim)
integer(kind=iwp) :: ibas, iOrb, iStart, Jorb, Korb
real(kind=wp) :: A, F
integer(kind=iwp), parameter :: Occ = 1, Vir = 0

#include "macros.fh"
unused_var(Center)

!write(u6,*) 'nDim,nDim',nDim,nDim
!write(u6,*) Center
do iOrb=1,nDim
  if ((iRestrict == 1) .and. (itype(iOrb) == Vir)) cycle
  F = Zero
  if (S(iOrb,iOrb) > Zero) F = One/sqrt(S(iOrb,iOrb))

  do ibas=1,nDim
    C(Ibas,Iorb) = F*C(Ibas,Iorb)
    !write(u6,*) C(Ibas,Iorb),IBas,Iorb
  end do

  do Jorb=1,nDim
    S(Iorb,Jorb) = F*S(Iorb,Jorb)
  end do
  do Jorb=1,nDim
    S(Jorb,Iorb) = F*S(Jorb,Iorb)
  end do

  iStart = 1
  if (iRestrict == 0) iStart = iOrb+1
  do jOrb=iStart,nDim
    if ((iRestrict == 1) .and. (itype(jOrb) == Occ)) cycle
    !if (Center(iOrb) == Center(jOrb)) cycle
    A = S(Iorb,Jorb)
    do Ibas=1,nDim
      C(Ibas,Jorb) = C(Ibas,Jorb)-A*C(Ibas,Iorb)
    end do
    do Korb=1,nDim
      S(Jorb,Korb) = S(Jorb,Korb)-A*S(Iorb,Korb)
    end do
    do Korb=1,nDim
      S(Korb,Jorb) = S(Korb,Jorb)-A*S(Korb,Iorb)
    end do
  end do
end do

return

end subroutine GramSchmidt
