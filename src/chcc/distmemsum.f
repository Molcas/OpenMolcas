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
        subroutine DistMemSum (NvGrp,maxdim,
     c        PossV1,PossV2,PossV3,
     c        PossH1,PossH2,
     c        PossT)

c
c       This routine do:
c       define initial possitions of V and H
c       used in summary routine
c

        implicit none
#include "chcc1.fh"
c
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3
        integer PossH1,PossH2
        integer PossT
c
c
c       help variables
        integer length
c
c
c1        V files
c
        length=no*no*maxdim*maxdim
        PossV1=possT
        PossT=PossT+length
        PossV2=possT
        PossT=PossT+length
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,99) 'DM V ',PossV1,PossV2,PossV3,length
        end if
c
c
c2        H files
c
        length=no*maxdim
        PossH1=possT
        PossT=PossT+length
        PossH2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,99) 'DM H ',PossH1,PossH2,PossV3,length
        end if
c
c
99      format (a7,10(i10,1x))
c
        if (printkey.ge.10) then
        write (6,99) 'PossT ',PossT
        end if
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(NvGrp)
        end
