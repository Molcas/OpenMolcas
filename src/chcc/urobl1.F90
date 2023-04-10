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

subroutine UrobL1(L1,NaGrp,LunAux)
! vyraba fily so simulovanymi L1 vektormi
! so spravnou strukturou

implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
integer NaGrp, LunAux
real*8 L1(1)
! help variables
integer aGrp, len_
real*8 schem

!1 cycle over a,be Groups

do aGrp=1,NaGrp

  !1.1 def length
  len_ = nc*DimGrpv(aGrp)*no

  !1.2 full L1 with random numbers
  schem = 1.0d-2
  call RNFill(len_,L1(1),schem)

  !1.3 open proper file
  !open(unit=LunAux,file=L1Name(aGrp),form='unformatted')
  call MOLCAS_BinaryOpen_Vanilla(LunAux,L1Name(aGrp))

  !1.4 write L1 into proper file
  write(6,*) aGrp,len_
  call wri_chcc(LunAux,len_,L1(1))

  close(LunAux)

end do

return

end subroutine UrobL1
