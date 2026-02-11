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

use Molcas, only: LenIn, MxAtom, MxSym
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MxBasis = 5000
integer(kind=iwp) :: iPrFmt, nBas(MxSym), nDel(MxSym), nNuc, nOcc(MxSym), nSym, nVir(MxSym)
real(kind=wp) :: GapThr, PrThr, SThr, TThr
character(len=LenIn+8) :: Label(MxBasis)
character(len=LenIn) :: AtName(MxAtom)
logical(kind=iwp) :: PrintEor, PrintMOs, PrintPop
#ifdef _HDF5_
integer(kind=iwp) :: wfn_energy, wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif
#ifdef _OLD_
real(kind=wp) :: xCharge(MxAtom)
#endif

public :: AtName, GapThr, iPrFmt, Label, MxBasis, nBas, nDel, nNuc, nOcc, nSym, nVir, PrintEor, PrintMOs, PrintPop, PrThr, SThr, &
          TThr
#ifdef _HDF5_
public :: wfn_energy, wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif
#ifdef _OLD_
public :: xCharge
#endif

end module GuessOrb_Global
