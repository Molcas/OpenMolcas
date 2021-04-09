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
        subroutine MakeWw (Ww,W1,W2,aSGrp,bSGrp,beSGrp,gaSGrp,key)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c                                     and W2(b",be",a",ga")
c        where index ab" is always a">b" (not a"=b")
c        and   index bega" is be">=ga" for W+ and be">ga" for W-
c       for all cases aSGrp >= bSGrp
c       N.B. for a"=b" W2 is not used (although it is passed through
c            the header (because of uniform shape) but no matter what
c            is there)
c       N.B. rutinky MakeWwHlpx niesu prilis vymakane @
c
c       parameter description:
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       W2     - array for W1(b",be",a",ga") (I)
c       xSGrp  - SubGroup of a",b",be",ga" (I)
c       key    - 1 - calc Ww+, 2 - calc Ww- (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 Ww(1)
        real*8 W1(1)
        real*8 W2(1)
        integer aSGrp,bSGrp,beSGrp,gaSGrp,key
c
c       help variables
        integer dimbe,dimga,dimbega,dima,dimb,dimab
c
c
c1      def dimensions
c
        dimbe=DimSGrpbe(beSGrp)
        dimga=DimSGrpbe(gaSGrp)
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpa(bSGrp)
c
        if (beSGrp.eq.gaSGrp) then
          if (key.eq.1) then
            dimbega=dimbe*(dimbe+1)/2
          else
            dimbega=dimbe*(dimbe-1)/2
          end if
        else
          dimbega=dimbe*dimga
        end if
c
        if (aSGrp.eq.bSGrp) then
          dimab=dima*(dima-1)/2
        else
          dimab=dima*dimb
        end if
c
c
c2      Make Ww matrix
c
        if (aSGrp.eq.bSGrp) then
          if (beSGrp.eq.gaSGrp) then
            call MakeWwHlp1 (Ww(1),W1(1),
     c                       dima,dimb,dimab,dimbe,dimga,dimbega,key)
          else
            call MakeWwHlp2 (Ww(1),W1(1),
     c                       dima,dimb,dimab,dimbe,dimga,key)
          end if
        else
          if (beSGrp.eq.gaSGrp) then
            call MakeWwHlp3 (Ww(1),W1(1),W2(1),
     c                       dima,dimb,dimbe,dimga,dimbega,key)
          else
            call MakeWwHlp4 (Ww(1),W1(1),W2(1),
     c                       dima,dimb,dimbe,dimga,key)
          end if
        end if
c
c
        return
        end
