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

subroutine SroMm1(nHer,MmSroG,la,lb,lr)
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

use Basis_Info, only: dbsc, nCnttp, Shells
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nHer, MmSroG, la, lb, lr
integer(kind=iwp) :: iacore, iAng, iCnttp, icoreb, icores, ip, iShll, ld, nac, ncb, nExpi, nOrder
! Statement function
integer(kind=iwp) :: nElem, i
nElem(i) = (i+1)*(i+2)/2

nOrder = 0
ld = 1
MmSroG = 0
do iCnttp=1,nCnttp
  if (.not. dbsc(iCnttp)%ECP) cycle
  do iAng=0,dbsc(iCnttp)%nSRO-1
    iShll = dbsc(iCnttp)%iSRO+iAng
    nExpi = Shells(iShll)%nExp
    if (nExpi == 0) cycle

    ip = 0

    nac = nElem(la)*nElem(iAng)
    ncb = nElem(iAng)*nElem(lb)

    ip = ip+6*nelem(la)*nelem(lb) ! final
    ip = ip+4*nac*nExpi ! FA1
    ip = ip+4*ncb*nExpi !FB1
    ip = ip+nExpi*nExpi !Tmp core
    ip = ip+nExpi !Tmp in sro

    nHer = (la+1+iAng+1+ld)/2
    nOrder = max(nHer,nOrder)
    iacore = 6+3*nHer*(la+1+ld)+3*nHer*(iAng+1)+3*nHer*(lr+1)+3*(la+1+ld)*(iAng+1)*(lr+1)+1

    nHer = (lb+1+iAng+1+ld)/2
    nOrder = max(nHer,nOrder)
    icoreb = 6+3*nHer*(lb+1+ld)+3*nHer*(iAng+1)+3*nHer*(lr+1)+3*(lb+1+ld)*(iAng+1)*(lr+1)+1

    icores = max(icoreb,iacore)*nExpi

    MmSroG = max(MmSroG,ip+icores)

  end do
end do
nHer = nOrder

return

end subroutine SroMm1
