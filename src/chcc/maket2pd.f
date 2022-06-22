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
        subroutine MakeT2pd (T2p,Tau,aGrp,aSGrp)
c
c       this routine do:
c       Make T2p((a>b)",i>=j )  = Tau(i,j,(a>=b)")+Tau(j,i,a>=b)")
c                           from  Tau((a>=b)',i,j)
c
c       parameter description:
c       T2p    - T2+ array (O)
c       Tau    - Tau array (I)
c       xGrp   - Group of a (I)
c       xSGrp  - SubGroup of a (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 T2p(1)
        real*8 Tau(1)
        integer aGrp,aSGrp
c
c       help variables
        integer dimi,dimij,dimap,dimapp,dimabp
c
        dimi=no
        dimij=no*(no+1)/2
c
        dimap=DimGrpa(aGrp)
        dimabp=dimap*(dimap+1)/2
c
        dimapp=DimSGrpa(aSGrp)
c
c
        call makeT2pdHlp (T2p(1),Tau(1),aGrp,aSGrp,
     c                   dimi,dimij,dimapp,dimap,dimabp)
c
        return
        end
