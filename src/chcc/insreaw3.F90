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
        subroutine InsReaW3 (aSGrp,beSGrp,bSGrp,length)
!
!        this routine do:
!        Check which W3 file corresponds to given combination
!        of indexes,
!        - increase the parameter, checking the overal number of
!        required integrals if these block is conted first time
!        - and set corresponding InqW3 to T
!
!
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "parcc.fh"
!
        integer aSGrp,beSGrp,bSGrp,length
!
!        help variables
!
        integer dima,dimb,dimbe
        integer abeSGrp
!
!
!0        def dimensions
!
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpa(bSGrp)
        dimbe=DimSGrpbe(beSGrp)
!
        if (aSGrp.gt.beSGrp) then
!1        case aSGrp>beSGrp, integrals (a",be"|b",i) to be read
          abeSGrp=(aSGrp*(aSGrp-1)/2)+beSGrp
          if (InqW3(abeSGrp,bSGrp).eqv..False.) then
            InqW3(abeSGrp,bSGrp)=.True.
            length=length+dima*dimbe*dimb*no
          end if
!
        else if (aSGrp.eq.beSGrp) then
!2        case aSGrp=beSGrp, integrals (a">=be"|b",i)to be read and expand
          abeSGrp=(aSGrp*(aSGrp-1)/2)+beSGrp
          if (InqW3(abeSGrp,bSGrp).eqv..False.) then
            InqW3(abeSGrp,bSGrp)=.True.
            length=length+(no*dima*(dima+1)*dimb)/2
          end if
!
        else
!3        case aSGrp<beSGrp, integrals (be",a"|b",i) to be read and mapped
          abeSGrp=(beSGrp*(beSGrp-1)/2)+aSGrp
          if (InqW3(abeSGrp,bSGrp).eqv..False.) then
            InqW3(abeSGrp,bSGrp)=.True.
            length=length+dima*dimbe*dimb*no
          end if
        end if
!
        return
        end
