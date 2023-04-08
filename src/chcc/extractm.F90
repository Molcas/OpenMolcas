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
        subroutine ExtractM (M,L2,cGrp,deGrp,cSGrp,deSGrp)
!
!       this routine do:
!       extract M(m,c",de") from L2(m,c',de')
!
!       parameter description:
!       M     - M file (O)
!       L     - L2 file (O)
!       xGrp  - c,delta Group (I)
!       xSGrp - c,delta Group (I)
!
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
!
        real*8 M(*),L2(*)
        integer cGrp,deGrp,cSGrp,deSGrp
!
!       help variables
        integer lenMCpp
        integer possM0,possL20,incL20
        integer i,k,depp
!
!
!1      Initial settings
!
!1.1    def length of mc" block (=increment for possM0)
        lenMCpp=nc*DimSGrpa(cSGrp)
!
!1.2    def possM0
        possM0=1
!
!1.3    def possL20
        PossL20=1
        k=0
        if (deSGrp.gt.Grpbelow(deGrp)) then
          do i=Grpbelow(deGrp),deSGrp-1
            k=k+DimSGrpbe(i)
          end do
        end if
        possL20=PossL20+nc*DimGrpa(cGrp)*k
!
        k=0
        if (cSGrp.gt.Grpalow(cGrp)) then
          do i=Grpalow(cGrp),cSGrp-1
            k=k+DimSGrpa(i)
          end do
        end if
        possL20=PossL20+nc*k
!
!1.4    def increment for possL20
        incL20=dimGrpa(cGrp)*nc
!
!2      cycle over de"
        do depp=1,DimSgrpbe(deSGrp)
!
!2.1      copy block of #mc" size
          call mv0u (lenMCpp,L2(possL20),1,M(possM0),1)
!
!2.2      upgrade possM0 and possL20 for next use
          possM0=possM0+lenMCpp
          possL20=possL20+incL20
!
        end do
!
!
        return
        end
