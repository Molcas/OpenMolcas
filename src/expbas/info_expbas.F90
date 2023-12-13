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

module info_expbas_mod

use Definitions, only: iwp

implicit none
private

#include "Molcas.fh"
! Exporting some useful parameters
public :: LenIn, MxAtom, mxsym

!> Different kinds of orbitals
!> f, i, 1, 2, 3, s, d
integer(kind=iwp), parameter :: n_orb_kinds = 7

integer(kind=iwp) :: nBas1(MxSym), nBas2(MxSym), nSym1, nSym2
#ifdef _HDF5_
integer(kind=iwp) :: wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif
logical(kind=iwp) :: DoExpbas, DoDesy
character(len=512) :: EB_FileOrb

public :: DoDesy, DoExpbas, EB_FileOrb, n_orb_kinds, nBas1, nBas2, nSym1, nSym2
#ifdef _HDF5_
public :: wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif

end module info_expbas_mod
