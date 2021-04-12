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
        subroutine ReaW (W,aSGrp,beSGrp,bSGrp,gaSGrp,LunInt)
c
c       simuluje citanie VVVV integralov
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 W(1)
        integer aSGrp,beSGrp,bSGrp,gaSGrp,LunInt
c
c       help variables
        integer dim
c
        dim=DimSGrpa(aSGrp)*DimSGrpbe(beSGrp)*
     c      DimSGrpa(bSGrp)*DimSGrpbe(gaSGrp)
c
        call rea1 (LunInt,dim,W(1))
c
        return
        end
