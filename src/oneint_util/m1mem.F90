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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

subroutine M1Mem( &
#                define _CALLING_
#                include "mem_interface.fh"
                )
!***********************************************************************
!  Object: to compute the number of reals the kernel routine will      *
!          need for the computation of a matrix element between two    *
!          cartesian Gaussian functions with the total angular momentum*
!          of la and lb (la=0 s-function, la=1 p-function, etc.)       *
!          lr is the order of the operator (this is only used when the *
!          integrals are computed with the Hermite-Gauss quadrature).  *
!                                                                      *
!  Called from: OneEl                                                  *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri3_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: iAng(4), k, MemM10, MemPrm, nFlop, nMem

#include "macros.fh"
unused_var(lr)

call mHrr(la,lb,nFlop,nMem)

iAng(1) = la
iAng(2) = lb
iAng(3) = 0
iAng(4) = 0
call MemRys(iAng,MemPrm)
MemM10 = 6+MemPrm
nHer = (la+lb+2)/2

k = nTri3_Elem1(la+lb)-nTri3_Elem1(max(la,lb)-1)

! nMem : memory for Hrr
! k    : memory for primitives in M1Int0
! MemM10: scratch in M1Int0

Mem = MemM10+max(nMem,k)

return

end subroutine M1Mem
