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
        subroutine Chck_B (BB,
     c             dima,dimb,dimbe,dimga,adda,addb,addbe,addga)
c
c        this routine test L2
c
        implicit none
#include "chcc1.fh"
        integer dima,dimb,dimbe,dimga,adda,addb,addbe,addga
         real*8 BB(1:dima,1:dimbe,1:dimb,1:dimga)
c
c        help var
        integer a,b,be,ga,bad,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do ga=addga+1,addga+dimga
        do b=addb+1,addb+dimb
        do be=addbe+1,addbe+dimbe
        do a=adda+1,adda+dima
c
          s=Bc(a,b,be,ga)
c
           if (abs(BB(a-adda,be-addbe,b-addb,ga-addga)-s)
     &              .gt.1.0d-10) then
            bad=bad+1
            BB(a-adda,be-addbe,b-addb,ga-addga)=s
          end if
          ntot=ntot+1

        end do
        end do
        end do
        end do
c
        write (6,*) ' B test ',bad,ntot
c
        return
        end
