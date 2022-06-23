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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine SROMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )
!***********************************************************************
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
#include "mem_interface.fh"
integer(kind=iwp) :: iAng, iCnttp, ip, iShll, memMlt, nac, ncb, nExpi, nH
! Statement function
integer(kind=iwp) :: nElem, i
nElem(i) = (i+1)*(i+2)/2

nHer = 0
Mem = 0
do iCnttp=1,nCnttp
  if (.not. dbsc(iCnttp)%ECP) cycle
  do iAng=0,dbsc(iCnttp)%nSRO-1
    iShll = dbsc(iCnttp)%iSRO+iAng
    nExpi = Shells(iShll)%nExp
    if (nExpi == 0) cycle

    ip = 0
    ip = ip+nExpi**2
    nac = nElem(la)*nElem(iAng)
    ip = ip+nExpi*nac
    ip = ip+3*nExpi
    ip = ip+nExpi
    ip = ip+nExpi
    ip = ip+nExpi

    call MltMmP(nH,MemMlt,la,iAng,lr)
    nHer = max(nH,nHer)
    Mem = max(Mem,ip+nExpi*MemMlt)
    ip = ip-6*nExpi

    ncb = nElem(iAng)*nElem(lb)
    ip = ip+nExpi*ncb
    ip = ip+3*nExpi
    ip = ip+nExpi
    ip = ip+nExpi
    ip = ip+nExpi

    call MltMmP(nH,MemMlt,iAng,lb,lr)
    nHer = max(nH,nHer)
    Mem = max(Mem,ip+nExpi*MemMlt)
    ip = ip-6*nExpi

    ip = ip+max(nExpi*nac,ncb*nExpi)
    Mem = max(Mem,ip)

  end do
end do

return

end subroutine SROMem
