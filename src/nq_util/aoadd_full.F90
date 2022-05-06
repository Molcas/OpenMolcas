!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991,2021,2022, Roland Lindh                           *
!***********************************************************************

subroutine AOAdd_Full(PrpInt,nPrp,nD)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January 1991                                             *
!***********************************************************************

use nq_Grid, only: iBfn_Index
use nq_Grid, only: AOInt => Dens_AO
implicit real*8(A-H,O-Z)
#include "real.fh"
real*8 PrpInt(nPrp,nD)
! Statement function
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *

nBfn = size(iBfn_index,2)
do iBfn=1,nBfn
  Indi = iBfn_Index(1,iBfn)

  do jBfn=1,iBfn
    Indj = iBfn_Index(1,jBfn)

    ! Add one matrix element

    PrpInt(iTri(Indi,Indj),:) = PrpInt(iTri(Indi,Indj),:)+AOInt(iBfn,jBfn,:)
  end do
end do

return

end subroutine AOAdd_Full
