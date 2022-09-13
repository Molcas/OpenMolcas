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

subroutine PrjMm1( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )
!***********************************************************************
!                                                                      *
!  Object: to compute the number of real*8 the kernel routine will     *
!          need for the computation of a matrix element between two    *
!          cartesian Gaussian functions with the total angular momentum*
!          of la and lb (la=0 s-function, la=1 p-function, etc.)       *
!          lr is the order of the operator (this is only used when the *
!          integrals are computed with the Hermite-Gauss quadrature).  *
!                                                                      *
!  Called from: OneEl                                                  *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: dbsc, nCnttp, Shells
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: iacore, iAng, iCnttp, icoreb, icores, ip, iShll, ld, nac, nBasisi, ncb, nExpi, nOrder

nOrder = 0
ld = 1
Mem = 0
do iCnttp=1,nCnttp
  if (.not. dbsc(iCnttp)%ECP) cycle
  do iAng=0,dbsc(iCnttp)%nPrj-1
    iShll = dbsc(iCnttp)%iPrj+iAng
    nExpi = Shells(iShll)%nExp
    nBasisi = Shells(iShll)%nBasis
    if ((nExpi == 0) .or. (nBasisi == 0)) cycle

    ip = 0

    nac = nTri_Elem1(la)*nTri_Elem1(iAng)
    ncb = nTri_Elem1(iAng)*nTri_Elem1(lb)

    ip = ip+6*nTri_Elem1(la)*nTri_Elem1(lb) ! final
    ip = ip+4*nac*nExpi ! FA1
    ip = ip+4*ncb*nExpi !FB1
    ip = ip+nExpi*nExpi !Tmp

    nHer = (la+1+iAng+1+ld)/2
    nOrder = max(nHer,nOrder)
    iacore = 6+3*nHer*(la+1+ld)+3*nHer*(iAng+1)+3*nHer*(lr+1)+3*(la+1+ld)*(iAng+1)*(lr+1)+1

    nHer = (lb+1+iAng+1+ld)/2
    nOrder = max(nHer,nOrder)
    icoreb = 6+3*nHer*(lb+1+ld)+3*nHer*(iAng+1)+3*nHer*(lr+1)+3*(lb+1+ld)*(iAng+1)*(lr+1)+1

    icores = max(icoreb,iacore)*nExpi

    Mem = max(Mem,ip+icores)

  end do
end do
nHer = nOrder
!write(u6,*) 'mem',Mem

return

end subroutine PrjMm1
