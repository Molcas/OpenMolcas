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
! Copyright (C) Yannick Carissan                                       *
!               2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GenerateP(Ovlp,cMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Debug)
! Author: Yannick Carissan.
!
! Modifications:
!    - October 6, 2005 (Thomas Bondo Pedersen):
!      Reduce operation count and use BLAS.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: nBasis, nOrb2Loc, nAtoms, nBas_per_Atom(*), nBas_Start(*)
real(kind=wp) :: Ovlp(nBasis,nBasis), cMO(nBasis,*), PA(nOrb2Loc,nOrb2Loc,nAtoms)
character(len=LenIn8) :: BName(*)
logical(kind=iwp) :: Debug
real(kind=wp), allocatable :: SBar(:,:)

call mma_Allocate(SBar,nBasis,nOrb2Loc,Label='SBar')
call GenerateP_1(Ovlp,cMO,Sbar,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Debug)
call mma_deallocate(SBar)

end subroutine GenerateP
