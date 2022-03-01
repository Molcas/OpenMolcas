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

subroutine zasun(i1,length,valn,jn,kn,ln)
! control routine over zasun process
!
! i1 - number of pivot index (I)
! length - number of valid integrals in block (I)

integer length, i1
#include "reorg.fh"
real*8 valn(1:nsize,1:mbas)
integer jn(1:nsize,1:mbas)
integer kn(1:nsize,1:mbas)
integer ln(1:nsize,1:mbas)

if (zrkey == 1) then
  call zasun_zr(i1,length,valn,jn,kn,ln)
else
  call zasun_pck(i1,length,valn,jn,kn,ln)
end if

return

end subroutine zasun
