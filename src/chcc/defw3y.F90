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
        subroutine  DefW3y (aGrp,bGrp,cGrp,w3y)
!
!        this routine do:
!        define w3y key, which indicates, if atleast one
!        W3 file need to be calculated on this node for given
!        a',b',c'
!
        implicit none
        integer aGrp,bGrp,cGrp,w3y
#include "chcc1.fh"
#include "o2v4.fh"
#ifdef _MOLCAS_MPP_
#include "parcc.fh"
!
!        help variables
        integer aSGrp,bSGrp,abSGrp,bSGrpUp,cSGrp
!
        w3y=0
!
!       cycle over a">=b" subgroups for a',b'
        do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
        if (aGrp.eq.bGrp) then
          bSGrpUp=aSGrp
        else
          bSGrpUp=GrpaUp(bGrp)
        end if
        do bSGrp=GrpaLow(bGrp),bSGrpUp
        abSGrp=aSGrp*(aSGrp-1)/2+bSGrp
!
!          cycle over c" for c'
          do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
          if (InqW3(abSGrp,cSGrp)) then
            w3y=w3y+1
          end if
          end do
!
!       end cycle over a">=b" subgroups
        end do
        end do
!
#else
        w3y=1
! Avoid unused argument warnings
        if (.false.) then
          call Unused_integer(aGrp)
          call Unused_integer(bGrp)
          call Unused_integer(cGrp)
        end if
#endif
!
        return
        end
