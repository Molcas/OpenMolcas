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

subroutine UrobTau(Tau,NaGrp,LunTau)
! vyraba subor LunTau so simulovanymi Tau amplitudami
! so spravnou strukturou

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: Tau(1)
integer(kind=iwp) :: NaGrp, LunTau
#include "chcc1.fh"
#include "o3v3.fh"
integer(kind=iwp) :: aGrp, bGrp, length

!1 cycle over Groups a,b

do aGrp=1,NaGrp
  do bGrp=1,aGrp

    !1.1 def legth

    if (aGrp == bGrp) then
      ! groups of a and b are equal, reading for a'>=b'
      length = no*no*DimGrpv(aGrp)*(DimGrpv(bGrp)+1)/2
    else
      ! aGrp>bGrp, reading for a',b' in given groups
      length = no*no*DimGrpv(aGrp)*DimGrpv(bGrp)
    end if

    !1.2 full Tau with random numbers

    call RNFill(length,Tau(1),1.0e-2_wp)

    !1.3 read block

    write(u6,*) aGrp,bGrp,length
    call wri_chcc(lunTau,length,Tau(1))

  end do
end do

!2 rewind file
rewind(LunTau)

return

end subroutine UrobTau
