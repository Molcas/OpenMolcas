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
        subroutine  DefW4y2 (aSGrp,bSGrp,cGrp,dGrp,w4y)
!
!        this routine do:
!        define w4y key, which indicates, if atleast one
!        W4 file need to be calculated on this node for
!        given a",b",c',d'
!
        implicit none
        integer aSGrp,bSGrp,cGrp,dGrp,w4y
#include "chcc1.fh"
#include "o2v4.fh"
#ifdef _MOLCAS_MPP_
#include "parcc.fh"
!
!        help variables
        integer abSGrp,cSGrp,cdSGrp,dSGrp,dSGrpUp
!
        w4y=0
!
!       def abSGrp
        abSGrp=aSGrp*(aSGrp-1)/2+bSGrp
!
!          cycle over  c">=d" for c',d'
          do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
          if (cGrp.eq.dGrp) then
            dSGrpUp=cSGrp
          else
            dSGrpUp=GrpaUp(dGrp)
          end if
          do dSGrp=GrpaLow(dGrp),dSGrpUp
          cdSGrp=cSGrp*(cSGrp-1)/2+dSGrp
            if (InqW4(abSGrp,cdSGrp)) then
              w4y=w4y+1
            end if
          end do
          end do
!
!
#else
        w4y=1
! Avoid unused argument warnings
        if (.false.) then
          call Unused_integer(aSGrp)
          call Unused_integer(bSGrp)
          call Unused_integer(cGrp)
          call Unused_integer(dGrp)
        end if
#endif
!
        return
        end
