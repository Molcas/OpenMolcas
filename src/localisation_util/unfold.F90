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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine UnFold(A,nAA,B,nBB,nSym,nBas)
!***********************************************************************
!                                                                      *
!     purpose: Expand the density matrix to full storage and scale     *
!              off-diagonal elements                                   *
!                                                                      *
!                                                                      *
!     input:                                                           *
!       A       : triangular matrix of length nAA                      *
!       nSym    : number of blocks to be expanded                      *
!       nBas    : leading dimensions in each block                     *
!                                                                      *
!     output:                                                          *
!       B       : expanded matrix of length nBB                        *
!                                                                      *
!     called from: PMat                                                *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 A(nAA), B(nBB)
integer nBas(nSym)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

k1 = 0
k2 = 0
Factor = Half
do iSym=1,nSym
  nb = nBas(iSym)
  do ib=1,nb
    do jb=1,(ib-1)
      l1 = jb+(ib-1)*nb+k2
      l2 = ib+(jb-1)*nb+k2
      l3 = jb+ib*(ib-1)/2+k1
      B(l1) = Factor*A(l3)
      B(l2) = Factor*A(l3)
    end do
    l2 = ib+(ib-1)*nb+k2
    l3 = ib+ib*(ib-1)/2+k1
    B(l2) = A(l3)
  end do
  k1 = k1+nb*(nb+1)/2
  k2 = k2+nb*nb
end do

end subroutine UnFold
