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

subroutine FMOVE_CVB(IA,IB,N)

implicit real*8(A-H,O-Z)
!integer N
!real*8 IA(N), IB(N)
real*8 IA, IB
dimension IA(*), IB(*)
#include "SysDef.fh"

!call DCOPY_(N*RtoI,IA,1,IB,1)
do I=1,N*RtoI
  IB(I) = IA(I)
end do

return

end subroutine FMOVE_CVB
