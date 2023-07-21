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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_PTS_Final(NVT,l_NVT)
!
! Thomas Bondo Pedersen, April 2010.

implicit none
integer l_NVT
integer NVT(l_NVT)
#include "choglob.fh"
integer i

call iCopy(l_NVT,NVT,1,NumCho_G,1)
NumChT_G = NumCho_G(1)
do i=2,l_NVT
  NumChT_G = NumChT_G+NumCho_G(i)
end do
call Cho_Final(.false.)

end subroutine Cho_PTS_Final
