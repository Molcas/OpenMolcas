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
        subroutine MakeWwd (Ww,W1,aSGrp,beSGrp,gaSGrp)
!
!       this routine do:
!       Make  Ww(+)((aa)",(bega)") from W1(a",be",a",ga")
!        and   index bega" is be">=ga"
!       for all cases aSGrp >= bSGrp
!       N.B. rutinky MakeWwHlpx niesu prilis vymakane @
!
!       parameter description:
!       Ww     - array for Ww+(-) (O)
!       W1     - array for W1(a",be",a",ga") (I)
!       xSGrp  - SubGroup of a",be",ga" (I)
!
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
!
        real*8 Ww(1)
        real*8 W1(1)
        integer aSGrp,beSGrp,gaSGrp
!
!       help variables
        integer dimbe,dimga,dimbega,dima
!
!
!1      def dimensions
!
        dimbe=DimSGrpbe(beSGrp)
        dimga=DimSGrpbe(gaSGrp)
        dima=DimSGrpa(aSGrp)
!
        if (beSGrp.eq.gaSGrp) then
          dimbega=dimbe*(dimbe+1)/2
        else
          dimbega=dimbe*dimga
        end if
!
!
!2      Make Ww matrix
!
        if (beSGrp.eq.gaSGrp) then
          call MakeWwdHlp1 (Ww(1),W1(1),                                &
     &                   dima,dimbe,dimbega)
        else
          call MakeWwdHlp2 (Ww(1),W1(1),                                &
     &                   dima,dimbe,dimga)
        end if

!
!
        return
        end
