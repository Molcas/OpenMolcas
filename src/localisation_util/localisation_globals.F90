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

use Molcas, only: MxSym, LenIn
use Definitions, only: wp, iwp
use Constants, only: Zero

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
! nOrb2Loc: Number of orbital to localise
!
! BName   : Basis function names
!
! LC_FileOrb: orbital file for Seward to read

integer(kind=iwp) :: AnalyseLoc, ChargeType, fileorb_id, inpOptMeth, iWave, LocModel, MoldMod, MxConstr, nActa, nAtoms, &
                     nBas(MxSym), nCMO, nConstr(MxSym), nFro(MxSym), nMxIter, nOccInp(MxSym), nOrb(MxSym), nOrb2Loc(MxSym), nSym, &
                     nVirInp(MxSym), OptMeth
#ifdef _HDF5_
integer(kind=iwp) :: wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif
real(kind=wp) :: bias = 10.0_wp, GEKThr_Grad = 0.1_wp, GEKThr_Kappa = 0.1_wp, ScrFac = Zero, ThrDomain(2), ThrGrad, &
                 ThrPairDomain(3), ThrRot, Thrs, ThrSel
logical(kind=iwp) :: AnaAtom, AnaDomain, Analysis, AnaPAO, AnaPAO_Save, ChoStart, Debug = .false., DoCNOs, DoDomain, EvalER, &
                     getIMmldn, isHDF5 = .false., LocCanOrb, LocNatOrb, LocPAO, Maximisation, Order, PrintMOs, Silent, Skip, &
                     Test_Localisation, Timing, useFH, Wave
character(len=512) :: LC_FileOrb
character(len=3) :: AnaNrm
integer(kind=iwp), allocatable :: Ind(:), nBas_per_Atom(:), nBas_Start(:), posel(:)
real(kind=wp), allocatable :: CMO(:), DispList(:,:), EOrb(:), FuncList(:), GradList(:,:), kappa_cnt(:,:), MOrig(:), Occ(:), &
                              Ovlp(:,:), Ovlp_sqrt(:,:), xkappa_cnt(:,:)
character(len=LenIn+8), allocatable :: BName(:)
character(len=LenIn), allocatable :: NamAct(:)

public :: AnaAtom, AnaDomain, AnalyseLoc, Analysis, AnaNrm, AnaPAO, AnaPAO_Save, bias, BName, ChargeType, ChoStart, CMO, Debug, &
          DispList, DoCNOs, DoDomain, EOrb, EvalER, fileorb_id, FuncList, GEKThr_Grad, GEKThr_Kappa, getIMmldn, GradList, Ind, &
          inpOptMeth, isHDF5, iWave, kappa_cnt, LC_FileOrb, LocCanOrb, LocModel, LocNatOrb, LocPAO, Maximisation, MoldMod, MOrig, &
          MxConstr, nActa, NamAct, nAtoms, nBas, nBas_per_Atom, nBas_Start, nCMO, nConstr, nFro, nMxIter, nOccInp, nOrb, nOrb2Loc, &
          nSym, nVirInp, Occ, OptMeth, Order, Ovlp, Ovlp_sqrt, posel, PrintMOs, ScrFac, Silent, Skip, Test_Localisation, &
          ThrDomain, ThrGrad, ThrPairDomain, ThrRot, Thrs, ThrSel, Timing, useFH, Wave, xkappa_cnt
#ifdef _HDF5_
public :: wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif

end module Localisation_globals
