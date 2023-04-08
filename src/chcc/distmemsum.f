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
        subroutine DistMemSum (NvGrp,maxdim,                            &
     &        PossV1,PossV2,PossV3,                                     &
     &        PossH1,PossH2,                                            &
     &        PossT)

!
!       This routine do:
!       define initial possitions of V and H
!       used in summary routine
!

        implicit none
#include "chcc1.fh"
!
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3
        integer PossH1,PossH2
        integer PossT
!
!
!       help variables
        integer length
!
!
!1        V files
!
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
!
!
!2        H files
!
        length=no*maxdim
        PossH1=possT
        PossT=PossT+length
        PossH2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,99) 'DM H ',PossH1,PossH2,PossV3,length
        end if
!
!
99      format (a7,10(i10,1x))
!
        if (printkey.ge.10) then
        write (6,99) 'PossT ',PossT
        end if
!
        return
! Avoid unused argument warnings
        if (.false.) Call Unused_integer(NvGrp)
        end
