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
! Copyright (C) 2015, Roland Lindh                                     *
!***********************************************************************

subroutine EMFInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the electromagnetic field radiation integrals     *
!         using a complex Gauss-Hermite quadrature.                    *
!                                                                      *
! Called from: OneEl                                                   *
!                                                                      *
! Calling    : RecPrt                                                  *
!              CCrtCmp                                                 *
!              CAssmbl                                                 *
!              CVelInt                                                 *
!              CCmbnVe                                                 *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemistry - Angstrom,             *
!             University of Uppsala, Sweden. December 2015             *
!***********************************************************************

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: iComp, iDCRT(0:7), ipA, ipAOff, ipAxyz, ipB, ipBOff, ipBxyz, ipQxyz, ipRes, iPrint, ipVxyz, iRout, &
                     iStabO(0:7), lDCRT, llOper, LmbdT, nDCRT, nip, nOp, nStabO
integer(kind=iwp), external :: NrOpr

#include "macros.fh"
unused_var(ZInv)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 195
iPrint = nPrint(iRout)

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+1+nOrdOp)*2
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+1+nOrdOp)*2
ipQxyz = nip
nip = nip+nZeta*3*(la+1+nOrdOp)*(lb+1+nOrdOp)*2
if (nOrdOp == 1) then
  ipVxyz = nip
  nip = nip+nZeta*6*(la+1)*(lb+1)*2
  ipA = nip
  nip = nip+nZeta
  ipB = nip
  nip = nip+nZeta
  ipRes = nip
  nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*nComp
else
  ipVxyz = nip
  ipA = nip
  ipB = nip
  ipRes = nip
  nip = nip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*nComp
end if
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'EMFInt: nip-1 > nArr*nZeta')
  write(u6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(u6,*) ' Abend in EMFInt'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In EMFInt: A',' ',A,1,3)
  call RecPrt(' In EMFInt: RB',' ',RB,1,3)
  call RecPrt(' In EMFInt: KVector',' ',CoorO,1,3)
  call RecPrt(' In EMFInt: P',' ',P,nZeta,3)
  write(u6,*) ' In EMFInt: la,lb=',la,lb
end if

rFinal(:,:,:,:) = Zero

call EMFInt_Internal(Array)

llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

do lDCRT=0,nDCRT-1

  ! Accumulate contributions

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array(ipRes),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,One)

end do

return

! This is to allow type punning without an explicit interface
contains

subroutine EMFInt_Internal(Array)

  real(kind=wp), target :: Array(*)
  complex(kind=wp), pointer :: zAxyz(:,:,:,:), zBxyz(:,:,:,:), zQxyz(:,:,:,:), zQxyz2(:,:,:,:,:), zVxyz(:,:,:,:,:)
  integer(kind=iwp) :: iBeta

  ! Compute the cartesian values of the basis functions angular part
  ! Note that these arrays are complex.

  call c_f_pointer(c_loc(Array(ipAxyz)),zAxyz,[nZeta,3,nHer,la+nOrdOp+1])
  call c_f_pointer(c_loc(Array(ipBxyz)),zBxyz,[nZeta,3,nHer,lb+nOrdOp+1])
  call c_f_pointer(c_loc(Array(ipQxyz)),zQxyz,[nZeta,3,la+nOrdOp+1,lb+nOrdOp+1])

  call CCrtCmp(Zeta,P,nZeta,A,zAxyz,la+nOrdOp,HerR(iHerR(nHer)),nHer,CoorO)
  call CCrtCmp(Zeta,P,nZeta,RB,zBxyz,lb+nOrdOp,HerR(iHerR(nHer)),nHer,CoorO)

  ! Compute the cartesian components for the multipole moment
  ! integrals. The integrals are factorized into components.

  call CAssmbl(zQxyz,zAxyz,la+nOrdOp,zBxyz,lb+nOrdOp,nZeta,HerW(iHerW(nHer)),nHer)

  ! Compute the cartesian components for the velocity integrals.
  ! The velocity components are linear combinations of overlap components.

  if (nOrdOp == 1) then
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

    call c_f_pointer(c_loc(Array(ipVxyz)),zVxyz,[nZeta,3,la+1,lb+1,2])
    call CVelInt(zVxyz,zQxyz,la,lb,Array(ipA),Array(ipB),nZeta)

    ! Combine the cartesian components to the full one electron integral.

    call CCmbnVe(zQxyz,nZeta,la,lb,Zeta,rKappa,Array(ipRes),nComp,zVxyz,CoorO,P)
    nullify(zVxyz)
  else
    call c_f_pointer(c_loc(Array(ipQxyz)),zQxyz2,[nZeta,3,la+1,lb+1,nOrdOp+1])
    call CCmbnMP(zQxyz2,nZeta,la,lb,nOrdOp,Zeta,rKappa,Array(ipRes),nComp,CoorO,P)
    nullify(zQxyz2)
  end if

  nullify(zAxyz,zBxyz,zQxyz)

  return

end subroutine EMFInt_internal

end subroutine EMFInt
