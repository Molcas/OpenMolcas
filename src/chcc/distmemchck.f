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
        subroutine  DistMemChck (PossV1,PossV2,PossV3,PossT)
c
c        Dimensions for DistMemCheck
c        dimensions: V1    - no*nv*nv2, nv2*nv2
c                    V2    - nc*nv2
c                    V3    - nc*nv*no

        implicit none
#include "chcc1.fh"
        integer PossV1,PossV2,PossV3,PossT
c
c        help variables
        integer len
c
        PossT=PossFree
c
        PossV1=PossT
        len=no*nv*nv*(nv+1)/2
        if ((nv*(nv+1)*nv*(nv+1)/4).gt.len) then
        len=nv*(nv+1)*nv*(nv+1)/4
        end if
        PossT=PossT+len
c
        PossV2=PossT
        len=nc*nv*(nv+1)/2
        PossT=PossT+len
c
        PossV3=PossT
        len=nc*no*nv
c
        PossT=PossT+len
c
        write (6,*) ' Poss ChCk',PossV1,PossV2,PossV3,PossT
c
        return
        end
