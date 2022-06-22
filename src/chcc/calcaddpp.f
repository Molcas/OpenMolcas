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
        subroutine CalcAddpp (aSGrp,addapp)
c
c        this routine calc addapp
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
        integer aSGrp,addapp
c
c        help var
        integer i
c
          addapp=0
        if (aSGrp.gt.1) then
          do i=1,aSGrp-1
            addapp=addapp+DimSGrpa(i)
          end do
        end if
c
        return
        end
