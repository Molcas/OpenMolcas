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
use Molcas, only: LenIn
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lthBas, nBas_Tot

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
nBas_tot = sum(nBas(1:nSym))
call mma_allocate(BName,nBas_tot,Label='BName')
call Get_cArray('Unique Basis Names',BName,(LenIn+8)*nBas_tot)
! read number of atoms
call Get_iScalar('Unique atoms',nAtoms)
! read nuclear potential
!call get_dScalar('PotNuc',PotNuc)
call Peek_dScalar('PotNuc',PotNuc)

! Compute lengths of matrices
lthBas = sum(nBas(1:nSym))
call mma_allocate(Atom,lthBas,Label='Atom')
call mma_allocate(BType,lthBas,Label='BType')

! Define atom and type
Atom(:) = BName(1:lthBas)(1:LenIn)
BType(:) = BName(1:lthBas)(LenIn+1:LenIn+8)

return

end subroutine R1IBas
