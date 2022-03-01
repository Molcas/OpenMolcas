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

subroutine getpp_zr(lunpublic,pp,ipp,length)

#include "SysDef.fh"
integer lunpublic, length
real*8 pp(1:length)
integer ipp(1:length)

read(lunpublic) pp,ipp

return

end subroutine getpp_zr
