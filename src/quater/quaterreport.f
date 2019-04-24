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
      subroutine quaterreport(V1best,V2best,V1,V2)
      implicit none
#include "geoms.fh"
#include "options.fh"
      real*8 V1best(3),V2best(3),V1(3),V2(3)
      real*8 v1dot,v2dot
      real*8 ddot_

      Write(6,*) "Number of geometries generated : ",ngeoms

      if (rotate) then
        v1dot=ddot_(3,V1best,1,V1,1)
        v2dot=ddot_(3,V2best,1,V2,1)
c
        Call Add_Info("V1_dot_product",[v1dot],1,8)
        Call Add_Info("V2_dot_product",[v2dot],1,8)
c
        Write(6,*) "V1best.V1 = ",v1dot
        Write(6,*) "V2best.V2 = ",v2dot
      end if

      end
