************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
        subroutine Calc_addG (aGrp,adda)
c
c        calc add constant from group
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
        integer aGrp,adda
c
c        help var
        integer i
c
        adda=0
        do i=1,agrp-1
          adda=adda+DimGrpa(i)
        end do
c
        return
        end
