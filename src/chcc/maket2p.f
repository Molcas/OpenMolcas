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
        subroutine MakeT2p (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,keyT)
c
c       this routine do:
c       Make T2p(i>=j,(a>b)" )  = Tau(i,j,(a>=b)")+Tau(j,i,a>=b)")
c                           from  Tau((a>=b)',i,j)
c        or Transposed (T(ab",ij)
c
c       parameter description:
c       T2p    - T2+ array (O)
c       Tau    - Tau array (I)
c       xGrp   - Group of a,b (I)
c       xSGrp  - SubGroup of a,b (I)
c        keyT   - 0 - make T(ij,ab")
c                1 - make T(ab",ij)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 T2p(1)
        real*8 Tau(1)
        integer aGrp,bGrp,aSGrp,bSGrp,keyT
c
c       help variables
        integer dimi,dimij,dimap,dimbp,dimapp,dimbpp,dimabp,dimabpp
c
        dimi=no
        dimij=no*(no+1)/2
c
        dimap=DimGrpa(aGrp)
        dimbp=DimGrpa(bGrp)
        if (aGrp.eq.bGrp) then
          dimabp=dimap*(dimap+1)/2
        else
          dimabp=dimap*dimbp
        end if
c
        dimapp=DimSGrpa(aSGrp)
        dimbpp=DimSGrpa(bSGrp)
        if (aSGrp.eq.bSGrp) then
          dimabpp=dimapp*(dimapp-1)/2
        else
          dimabpp=dimapp*dimbpp
        end if
c
        if (keyT.eq.0) then
c        T+(ij,ab") case
c
        if (aGrp.eq.bGrp) then
          if (aSGrp.eq.bSGrp) then
            call makeT2pHlp1 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                        dimi,dimij,dimapp,dimabpp,dimap,dimabp)
          else
            call makeT2pHlp2 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                        dimi,dimij,dimapp,dimbpp,dimapp,dimabp)
          end if
        else
          call makeT2pHlp3 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                      dimi,dimij,dimapp,dimbpp,dimap,dimbp)
        end if
c
        else
c        T+(ab",ij) case
c
        if (aGrp.eq.bGrp) then
          if (aSGrp.eq.bSGrp) then
            call makeT2ptHlp1 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                        dimi,dimij,dimapp,dimabpp,dimap,dimabp)
          else
            call makeT2ptHlp2 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                        dimi,dimij,dimapp,dimbpp,dimapp,dimabp)
          end if
        else
          call makeT2ptHlp3 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                      dimi,dimij,dimapp,dimbpp,dimap,dimbp)
        end if
c
        end if

c
        return
        end
