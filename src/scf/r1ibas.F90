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

subroutine R1IBas()
!***********************************************************************
!                                                                      *
!     purpose: Read basis set informations.                            *
!                                                                      *
!***********************************************************************

use InfSCF, only: Atom, BName, BType, Header, nAtoms, nBas, nSym, PotNuc
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: i, iSym, LthBas, nBas_Tot

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! read file header
call Get_cArray('Seward Title',Header,144)
! read number of symm. species
call Get_iScalar('nSym',nSym)
! read number of basis functions per symmetry species
call Get_iArray('nBas',nBas,nSym)
! read basis function labels
nBas_tot = 0
do iSym=1,nSym
  nBas_tot = nBas_tot+nBas(iSym)
end do
call mma_allocate(BName,nBas_tot,Label='BName')
call Get_cArray('Unique Basis Names',BName,(LenIn8)*nBas_tot)
! read number of atoms
call Get_iScalar('Unique atoms',nAtoms)
! read nuclear potential
!call get_dScalar('PotNuc',PotNuc)
call Peek_dScalar('PotNuc',PotNuc)

! Compute lengths of matrices
lthBas = 0
do iSym=1,nSym
  lthBas = lthBas+nBas(iSym)
end do
call mma_allocate(Atom,lthBas,Label='Atom')
call mma_allocate(BType,lthBas,Label='BType')

! Define atom and type
do i=1,lthBas
  Atom(i) = BName(i)(1:LenIn)
  BType(i) = BName(i)(LenIn1:LenIn8)
end do

return

end subroutine R1IBas
