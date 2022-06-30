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

subroutine PPMmG( &
#                define _CALLING_
#                include "mem_interface.fh"
                )

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: lalbm, lalbp, lambl, lapbl

#include "macros.fh"
unused_var(lr)

nHer = 0
Mem = 0

lapbl = max(nTri_Elem1(la+1),nTri_Elem1(lb))**2
Mem = Mem+2*lapbl

if (la > 0) then
  lambl = max(nTri_Elem1(la-1),nTri_Elem1(lb))**2
else
  lambl = 0
end if
Mem = Mem+2*lambl

lalbp = max(nTri_Elem1(la),nTri_Elem1(lb+1))**2
Mem = Mem+2*lalbp

if (lb > 0) then
  lalbm = max(nTri_Elem1(la),nTri_Elem1(lb-1))**2
else
  lalbm = 0
end if
Mem = Mem+2*lalbm

return

end subroutine PPMmG
