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
        subroutine MakeWwd (Ww,W1,aSGrp,beSGrp,gaSGrp)
c
c       this routine do:
c       Make  Ww(+)((aa)",(bega)") from W1(a",be",a",ga")
c        and   index bega" is be">=ga"
c       for all cases aSGrp >= bSGrp
c       N.B. rutinky MakeWwHlpx niesu prilis vymakane @
c
c       parameter description:
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",a",ga") (I)
c       xSGrp  - SubGroup of a",be",ga" (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 Ww(1)
        real*8 W1(1)
        integer aSGrp,beSGrp,gaSGrp
c
c       help variables
        integer dimbe,dimga,dimbega,dima
c
c
c1      def dimensions
c
        dimbe=DimSGrpbe(beSGrp)
        dimga=DimSGrpbe(gaSGrp)
        dima=DimSGrpa(aSGrp)
c
        if (beSGrp.eq.gaSGrp) then
          dimbega=dimbe*(dimbe+1)/2
        else
          dimbega=dimbe*dimga
        end if
c
c
c2      Make Ww matrix
c
        if (beSGrp.eq.gaSGrp) then
          call MakeWwdHlp1 (Ww(1),W1(1),
     c                   dima,dimbe,dimbega)
        else
          call MakeWwdHlp2 (Ww(1),W1(1),
     c                   dima,dimbe,dimga)
        end if

c
c
        return
        end
