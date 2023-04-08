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
        subroutine MakeT2pd (T2p,Tau,aGrp,aSGrp)
!
!       this routine do:
!       Make T2p((a>b)",i>=j )  = Tau(i,j,(a>=b)")+Tau(j,i,a>=b)")
!                           from  Tau((a>=b)',i,j)
!
!       parameter description:
!       T2p    - T2+ array (O)
!       Tau    - Tau array (I)
!       xGrp   - Group of a (I)
!       xSGrp  - SubGroup of a (I)
!
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
!
        real*8 T2p(1)
        real*8 Tau(1)
        integer aGrp,aSGrp
!
!       help variables
        integer dimi,dimij,dimap,dimapp,dimabp
!
        dimi=no
        dimij=no*(no+1)/2
!
        dimap=DimGrpa(aGrp)
        dimabp=dimap*(dimap+1)/2
!
        dimapp=DimSGrpa(aSGrp)
!
!
        call makeT2pdHlp (T2p(1),Tau(1),aGrp,aSGrp,                     &
     &                   dimi,dimij,dimapp,dimap,dimabp)
!
        return
        end
