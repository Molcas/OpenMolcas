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

subroutine UrobI1(I1,NaGrp,LunAux)
! vyraba fily so simulovanymi I1 vektormi
! so spravnou strukturou

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimGrpv, I1Name, no
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: I1(*)
integer(kind=iwp), intent(in) :: NaGrp, LunAux
integer(kind=iwp) :: aGrp, len_
real(kind=wp) :: schem

!1 cycle over a,be Groups

do aGrp=1,NaGrp

  !1.1 def length
  len_ = no*DimGrpv(aGrp)*nTri_Elem(no)

  !1.2 full I1 with random numbers
  schem = 1.0e-2_wp
  call RNFill(len_,I1,schem)

  !1.3 open proper file
  !open(unit=LunAux,file=I1Name(aGrp),form='unformatted')
  call MOLCAS_BinaryOpen_Vanilla(LunAux,I1Name(aGrp))

  !1.4 write I1 into proper file
  write(u6,*) aGrp,len_
  call wri_chcc(LunAux,len_,I1)

  close(LunAux)

end do

return

end subroutine UrobI1
