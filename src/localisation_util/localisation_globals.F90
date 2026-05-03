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
use Constants, only: Zero
use Molcas, only: MxSym, LenIn

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


! Variables for undershot avoidance
type Loosen_Type
  real(kind=wp) :: Thrs, Thrs2, Step, Factor
end type Loosen_Type


integer(kind=iwp) :: fileorb_id, iWave, LocModel, LuSpool, MxConstr, nActa, nAtoms, nBas(MxSym), nCMO, nConstr(MxSym), &
                     nFro(MxSym), nMxIter, nOccInp(MxSym), nOrb(MxSym), nOrb2Loc(MxSym), nSym, nVirInp(MxSym), OptMeth, &
                     ChargeType, LocOrb, AnalyseLoc, MoldMod, inpOptMeth
#ifdef _HDF5_
integer(kind=iwp) :: wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif
real(kind=wp) :: Thrs, ThrRot, ThrGrad, ThrDomain(2), ThrPairDomain(3), ThrSel, ThrStep
logical(kind=iwp) :: AnaAtom, AnaDomain, Analysis, AnaPAO, AnaPAO_Save, ChoStart, DoCNOs, DoDomain, EvalER, isHDF5 = .false., &
                     LocCanOrb, LocNatOrb, LocPAO, Maximisation, Order, PrintMOs, Silent, Skip, Test_Localisation, Timing,&
                     Wave,Thrs_UsrDef, LocModel_UsrDef, nFro_UsrDef, nOrb2Loc_UsrDef,Freeze, getIMmldn, Debug = .false., useFH
character(len=512) :: LC_FileOrb
character(len=3) :: AnaNrm
integer(kind=iwp), allocatable :: Ind(:),nBas_per_Atom(:), nBas_Start(:)
character(len=LenIn+8), allocatable :: BName(:)
character(len=LenIn), allocatable :: NamAct(:)
real(kind=wp), allocatable :: CMO(:), EOrb(:), MOrig(:), Occ(:), FuncList(:), GradList(:,:), DispList(:,:),UmatList(:,:,:),&
                              kappa_cnt(:,:), xkappa_cnt(:,:), Ovlp(:,:), Ovlp_sqrt(:,:)

real(kind=wp):: GEKThr_Kappa = 1.0_wp,&
                GEKThr_Grad = 10.0_wp,&
                bias = 100.0_wp,&
                SOFact = 1.0e9_wp

real(kind=wp) :: ScrFac=Zero
type(Loosen_Type) :: Loosen

public :: AnaAtom, AnaDomain, Analysis, AnaNrm, AnaPAO, AnaPAO_Save, BName, ChoStart, CMO, DoCNOs, DoDomain, EOrb, EvalER, &
          fileorb_id, Ind, isHDF5, iWave, LC_FileOrb, LocCanOrb, LocModel, LocNatOrb, LocPAO, LuSpool, Maximisation, MOrig, &
          MxConstr, nActa, NamAct, nAtoms, nBas, nCMO, nConstr, nFro, nMxIter, nOccInp, nOrb, nOrb2Loc, nSym, nVirInp, Occ, Order, &
          PrintMOs, Silent, Skip, Test_Localisation, ThrDomain, ThrGrad, ThrPairDomain, ThrRot, Thrs, ThrSel, Timing, Wave,&
          ScrFac,Debug, OptMeth, ChargeType, LocOrb, Thrs_UsrDef, LocModel_UsrDef, nFro_UsrDef, nOrb2Loc_UsrDef,&
          Freeze,Loosen,FuncList,GradList,DispList,UmatList, ThrStep, GEKThr_Kappa, GEKThr_Grad, bias, SOFact, AnalyseLoc,&
          kappa_cnt, xkappa_cnt, Ovlp, Ovlp_sqrt, nBas_per_Atom, nBas_Start, MoldMod, getIMmldn, useFH, inpOptMeth
#ifdef _HDF5_
public :: wfn_fileid, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
#endif

end module Localisation_globals
