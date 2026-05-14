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

subroutine UrobL2(L2,NaGrp,NbeGrp,LunAux)
! vyraba fily so simulovanymi Cholesky vektormi
! so spravnou strukturou

use chcc_global, only: DimGrpv, L2Name, nc
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: L2(*)
integer(kind=iwp), intent(in) :: NaGrp, NbeGrp, LunAux
integer(kind=iwp) :: aGrp, beGrp, len_
real(kind=wp) :: schem

!1 cycle over a,be Groups

do aGrp=1,NaGrp
  do beGrp=1,NbeGrp

    !1.1 def length
    len_ = nc*DimGrpv(aGrp)*DimGrpv(beGrp)

    !1.2 full L2 with random numbers
    schem = 1.0e-2_wp
    call RNFill(len_,L2,schem)

    !1.3 open proper file
    !open(unit=LunAux,file=L2Name(aGrp,beGrp),form='unformatted')
    call MOLCAS_BinaryOpen_Vanilla(LunAux,L2Name(aGrp,beGrp))

    !1.4 write L2 into proper file
    write(u6,*) aGrp,beGrp,len_
    call wri_chcc(LunAux,len_,L2)

    close(LunAux)

  end do
end do

return

end subroutine UrobL2
