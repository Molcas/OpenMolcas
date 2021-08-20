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

subroutine GramSchmidt(S,C,nDim,type,Center,iRestrict)

implicit real*8(A-H,O-Z)
#include "real.fh"
real*8 S(nDim,nDim), C(nDim,nDim)
integer type(nDim), Occ, Vir, Center(nDim)
parameter(Occ=1,Vir=0)

!write(6,*) 'nDim,nDim',nDim,nDim
!write(6,*) Center
do iOrb=1,nDim
  if ((iRestrict == 1) .and. (type(iOrb) == Vir)) Go To 99
  F = Zero
  if (S(iOrb,iOrb) > Zero) F = One/sqrt(S(iOrb,iOrb))

  do ibas=1,nDim
    C(Ibas,Iorb) = F*C(Ibas,Iorb)
    !write(6,*) C(Ibas,Iorb),IBas,Iorb
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
    if ((iRestrict == 1) .and. (type(jOrb) == Occ)) Go To 98
    !if (Center(iOrb) == Center(jOrb)) Go To 98
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
98  continue
  end do
99 continue
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(Center)

end subroutine GramSchmidt
