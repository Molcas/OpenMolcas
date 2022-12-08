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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine KnEInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the kinetic energy integrals with the Gauss-      *
!         Hermite quadrature.                                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Modified to multipole moments November '90               *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerW
use Index_Functions, only: nTri_Elem1
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "rmat_option.fh"
#include "rmat.fh"
#include "print.fh"
integer(kind=iwp) :: iBeta, icop, ipA, ipAOff, ipAxyz, ipB, ipBOff, ipBxyz, ipDi, ipqC, ipQxyz, iPrint, ipRnr, ipRxyz, ipTxyz, &
                     iRout, lsum, nip
logical(kind=iwp) :: ABeq(3)

#include "macros.fh"
unused_var(ZInv)
unused_var(lOper)
unused_var(iChO)
unused_var(iStabM)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 150
iPrint = nPrint(iRout)
ABeq(:) = A == RB

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+2)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+2)
ipRxyz = nip
nip = nip+nZeta*3*nHer*(nOrdOp-1)
ipQxyz = nip
nip = nip+nZeta*3*(la+2)*(lb+2)*(nOrdOp-1)
ipTxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+1)
ipA = nip
nip = nip+nZeta
ipB = nip
nip = nip+nZeta
!                                                                      *
!***********************************************************************
!                                                                      *
if (RMat_type_integrals) then
  ipRnr = nip
  nip = nip+nZeta*(la+lb+3)
  ipqC = nip
  nip = nip+nZeta*(la+lb+1)
  ipDi = nip
  nip = nip+nZeta*(la+lb+1)
else
  ipRnr = -1
  ipqC = -1
  ipDi = -1
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'KNEInt: nip-1 > nArr*nZeta')
  write(u6,*) 'nip=',nip
  write(u6,*) 'nArr,nZeta=',nArr,nZeta
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In KnEInt: A',' ',A,1,3)
  call RecPrt(' In KnEInt: RB',' ',RB,1,3)
  call RecPrt(' In KnEInt: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In KnEInt: P',' ',P,nZeta,3)
  write(u6,*) ' In KnEInt: la,lb=',la,lb
end if

if (RMat_type_integrals) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! R-matrix calculations: continuum basis functions (A=B=P=0)
  !  Compute the contributions of the basis functions and multipole
  !  radial part

  lsum = la+lb+2
  call radlc(Zeta,nZeta,lsum,Array(ipRnr))

  ! Optional for photoionization:
  ! R-matrix calculations: continuum basis functions (A=B=P=0)
  !  Compute the contributions of the Coulomb operator times qCoul
  !  outside the sphere Omega

  if (abs(qCoul) > Epsq) then
    lsum = la+lb
    icop = 1
    call radlq(Zeta,nzeta,lsum,Array(ipqC),icop)
  end if

  if (abs(dipol1) > Epsq) then
    lsum = la+lb
    icop = 2
    call radlq(Zeta,nzeta,lsum,Array(ipDi),icop)
  end if

  ! Combine the radial and angular component to the full one electron integral.

  call CmbnKEr(Array(ipRnr),Array(ipqC),Array(ipDi),nZeta,la,lb,Zeta,rFinal,nComp,Alpha,nAlpha,Beta,nBeta)

else
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  ! Compute the cartesian values of the basis functions angular part

  call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),la+1,HerR(iHerR(nHer)),nHer,ABeq)
  call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),nHer,ABeq)

  ! Compute the contribution from the multipole moment operator

  ABeq(:) = .false.
  call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),nOrdOp-2,HerR(iHerR(nHer)),nHer,ABeq)

  ! Compute the cartesian components for the multipole moment
  ! integrals. The integrals are factorized into components.

  call Assmbl(Array(ipQxyz),Array(ipAxyz),la+1,Array(ipRxyz),nOrdOp-2,Array(ipBxyz),lb+1,nZeta,HerW(iHerW(nHer)),nHer)

  ! Compute the cartesian components for the kinetic energy integrals.
  ! The kinetic energy components are linear combinations of overlap components.

  ipAOff = ipA-1
  do iBeta=1,nBeta
    Array(ipAOff+1:ipAOff+nAlpha) = Alpha
    ipAOff = ipAOff+nAlpha
  end do

  ipBOff = ipB-1
  do iBeta=1,nBeta
    Array(ipBOff+1:ipBOff+nAlpha) = Beta(iBeta)
    ipBOff = ipBOff+nAlpha
  end do

  call Kntc(Array(ipTxyz),Array(ipQxyz),la,lb,Array(ipA),Array(ipB),nZeta)

  ! Combine the cartesian components to the full one electron integral.

  call CmbnKE(Array(ipQxyz),nZeta,la,lb,nOrdOp-2,Zeta,rKappa,rFinal,nComp,Array(ipTxyz))

end if

return

end subroutine KnEInt
