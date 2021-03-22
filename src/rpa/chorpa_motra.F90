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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoRPA_MOTra(includeFrozen,includeDeleted)

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Transform Cholesky vectors to MO basis.
!
! TODO/FIXME:
! 1. This routine computes all MO blocks (ij, ai, ab), even though
!    we may only need some of them. A more flexible interface would
!    be nice to have. Presumably not a performance issue in RPA,
!    though (remains to be verified).
! 2. For unrestricted calculations, the alpha and beta
!    transformations are done separately, which means that the AO
!    vectors are read twice. Simultaneous transformation would be
!    desirable!

implicit none
logical includeFrozen, includeDeleted
#include "rpa_data.fh"
#include "WrkSpc.fh"

character*12 SecNam
parameter(SecNam='ChoRPA_MOTra')

integer RPA_iUHF
external RPA_iUHF

character*6 BName(2)

integer nSpin
integer iSym
integer ip_lCMO, l_lCMO
integer ip, l
!integer ip_lnBas
!integer ip_lnOrb
integer ip_lnFro
integer ip_lnOcc
integer ip_Zeros
integer ip_lnVir
integer ip_lnDel
integer iSpin

nSpin = RPA_iUHF()
if (nSpin == 1) then
  BName(1) = 'MOVECS'
  BName(2) = 'unused'
else if (nSpin == 2) then
  BName(1) = 'MOVECa'
  BName(2) = 'MOVECb'
else
  call RPA_Warn(3,SecNam//': illegal nSpin')
  BName(1) = 'unused'
  BName(2) = 'unused'
end if
l_lCMO = nBas(1)**2
do iSym=2,nSym
  l_lCMO = l_lCMO+nBas(iSym)**2
end do
call GetMem('locCMO','Allo','Real',ip_lCMO,l_lCMO)
l = 5*nSym
call GetMem('local','Allo','Inte',ip,l)
ip_lnFro = ip
ip_lnOcc = ip+nSym
ip_Zeros = ip+2*nSym
ip_lnVir = ip+3*nSym
ip_lnDel = ip+4*nSym
call iZero(iWork(ip_Zeros),nSym)
if (includeFrozen) ip_lnFro = ip_Zeros
if (includeDeleted) ip_lnDel = ip_Zeros

do iSpin=1,nSpin

  ! Set orbital blocks
  if (includeFrozen) then
    do iSym=1,nSym
      iWork(ip_lnOcc-1+iSym) = nFro(iSym,iSpin)+nOcc(iSym,iSpin)
    end do
  else
    call iCopy(nSym,nFro(1,iSpin),1,iWork(ip_lnFro),1)
    call iCopy(nSym,nOcc(1,iSpin),1,iWork(ip_lnOcc),1)
  end if
  if (includeDeleted) then
    do iSym=1,nSym
      iWork(ip_lnVir-1+iSym) = nVir(iSym,iSpin)+nDel(iSym,iSpin)
    end do
  else
    call iCopy(nSym,nVir(1,iSpin),1,iWork(ip_lnVir),1)
    call iCopy(nSym,nDel(1,iSpin),1,iWork(ip_lnDel),1)
  end if
  ! Reorder CMO array
  call ChoRPA_MOTra_ReorderCMO(nSym,nBas,nOrb,nFro(1,iSpin),nOcc(1,iSpin),nVir(1,iSpin),nDel(1,iSpin),Work(ip_CMO(iSpin)), &
                               Work(ip_lCMO))
  ! Set base name for MO files
  ! Transform Cholesky vectors
  call Cho_MOTra_Internal(Work(ip_lCMO),l_lCMO,nSym,nBas,nOrb,iWork(ip_lnFro),iWork(ip_lnOcc),iWork(ip_Zeros),iWork(ip_lnVir), &
                          iWork(ip_lnDel),BName(iSpin),.false.,0,.false.)

end do

call GetMem('local','Free','Inte',ip,l)
call GetMem('locCMO','Free','Real',ip_lCMO,l_lCMO)

end subroutine ChoRPA_MOTra
