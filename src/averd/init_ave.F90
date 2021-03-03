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

subroutine Init_ave(Title,iPrint,Wset,Wsum,PrOcc,PrEne,DensityBased,ThrOcc,Dummy,iDummy)

implicit real*8(a-h,o-z)

!-- Initializations and defaults.
#include "mxdm.fh"
#include "real.fh"
#include "mxave.fh"

dimension Wset(MxSets)
logical PrOcc, PrEne, DensityBased
character Title*72

iPrint = 2
do 10,i=1,MxSets
  Wset(i) = 0.0d0
10 continue
Wsum = 0.0d0
PrOcc = .true.
PrEne = .false.
DensityBased = .true.
ThrOcc = 1.0d-5
Dummy = 1.0
iDummy = 1
Title = ' Untitled job. '

return

end subroutine Init_ave
