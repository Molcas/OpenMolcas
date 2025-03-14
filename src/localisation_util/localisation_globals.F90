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

module Localisation_globals

use Definitions, only: wp, iwp

implicit none
private

! CMO     : the MO coefficients
! nCMO    : size of CMO
! Occ     : occupancy vector
! Eor     : orbital energy vector
! ipInd   : orbital type info vector
! nSym    : number of Irrep
! nBas    : number of basis function per Irrep
! nOrb    : number of occ orb per irrep
! LuSpool : Unit number of the input
! nOrb2Loc: Number of orbital to localise
!
! BName   : Basis function names
!
! LC_FileOrb: orbital file for Seward to read

#include "Molcas.fh"

integer(kind=iwp) :: fileorb_id, iWave, LocModel, LuSpool, MxConstr, nActa, nAtoms, nBas(MxSym), nCMO, nConstr(MxSym), &
                     nFro(MxSym), nMxIter, nOccInp(MxSym), nOrb(MxSym), nOrb2Loc(MxSym), nSym, nVirInp(MxSym)
#ifdef _HDF5_
integer(kind=iwp) :: wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif
real(kind=wp) :: Thrs, ThrRot, ThrGrad, ThrDomain(2), ThrPairDomain(3), ThrSel
logical(kind=iwp) :: AnaAtom, AnaDomain, Analysis, AnaPAO, AnaPAO_Save, ChoStart, DoCNOs, DoDomain, EvalER, isHDF5 = .false., &
                     LocCanOrb, LocNatOrb, LocPAO, Maximisation, Order, PrintMOs, Silent, Skip, Test_Localisation, Timing, Wave
character(len=512) :: LC_FileOrb
character(len=3) :: AnaNrm
integer(kind=iwp), allocatable :: Ind(:)
character(len=LenIn8), allocatable :: BName(:)
character(len=LenIn), allocatable :: NamAct(:)
real(kind=wp), allocatable :: CMO(:), EOrb(:), MOrig(:), Occ(:)

public :: AnaAtom, AnaDomain, Analysis, AnaNrm, AnaPAO, AnaPAO_Save, BName, ChoStart, CMO, DoCNOs, DoDomain, EOrb, EvalER, &
          fileorb_id, Ind, isHDF5, iWave, LC_FileOrb, LocCanOrb, LocModel, LocNatOrb, LocPAO, LuSpool, Maximisation, MOrig, &
          MxConstr, nActa, NamAct, nAtoms, nBas, nCMO, nConstr, nFro, nMxIter, nOccInp, nOrb, nOrb2Loc, nSym, nVirInp, Occ, Order, &
          PrintMOs, Silent, Skip, Test_Localisation, ThrDomain, ThrGrad, ThrPairDomain, ThrRot, Thrs, ThrSel, Timing, Wave
#ifdef _HDF5_
public :: wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif

end module Localisation_globals

