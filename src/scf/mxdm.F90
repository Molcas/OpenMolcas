!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!----------------------------------------------------------------------*
!     Define global parameters.                                        *
!----------------------------------------------------------------------*
!     MxSym  - maximal number of symmetries (D2h group)                *
!     MxBas  - maximal number of basis functions per symmetry          *
!     MxBS   - MxBas*MxSym                                             *
!     MxAtms - maximal number of symmetry unique atoms                 *
!     MxOptm - maximal number of density matrices we want to perform   *
!              optimization on                                         *
!     MxKeep - size of the Map vector ( = MxOptm)                      *
!     MxDDsk - maximal number of density matrices stored on the disk   *
!     MxTit  - maximal number of title lines                           *
!     MxKp2U - maximal number of terms in U=exp(K) expansion           *
!----------------------------------------------------------------------*

module MxDM

use Definitions, only: iwp

implicit none
private

#include "Molcas.fh"
integer(kind=iwp), parameter :: MxAtms = MxAtom, MxBS = MxBas*MxSym, MxIter = 400, MxKeep = MxIter, MxDDsk = MxKeep, MxKp2U = 100, &
                                MxOptm = 20, MxTit = 1

public :: LenIn, LenIn8, MaxBfn, MxDDsk, MxIter, MxKeep, MxOptm, MxOrb, MxRoot, MxSym, MxTit

end module MxDM
