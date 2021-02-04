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

module GuessOrb_Global

use Definitions, only: wp, iwp

implicit none
private

#include "Molcas.fh"
integer(kind=iwp), parameter :: MxBasis = 5000
real(kind=wp) :: xCharge(MxAtom), PrThr, SThr, TThr, GapThr
integer(kind=iwp) :: nSym, nBas(MxSym), nOcc(MxSym), nVir(MxSym), nDel(MxSym), nNuc, iPrFmt
character(len=LenIn) :: Name(MxAtom)
character(len=LenIn8) :: Label(MxBasis)
logical(kind=iwp) :: PrintMOs, PrintEor, PrintPop
#ifdef _HDF5_
integer(kind=iwp) :: wfn_fileid, wfn_energy, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif

public :: GapThr, iPrFmt, Label, LenIn, LenIn1, LenIn8, MxAtom, MxBasis, MxSym, Name, nBas, nDel, nNuc, nOcc, nSym, &
          nVir, PrintEor, PrintMOs, PrintPop, PrThr, SThr, TThr, xCharge
#ifdef _HDF5_
public :: wfn_energy, wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif

end module GuessOrb_Global
