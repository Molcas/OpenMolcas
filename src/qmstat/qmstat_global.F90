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

module qmstat_global

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

!----------------------------------------------------------------------*
! Common variables shared by all QmTypes.
!----------------------------------------------------------------------*
! INTEGER:
! --------
! lMax          -        How many bases in solvent region.
! info_atom     -        Atomic number of QM atoms.
!
! REAL:
! -----
! QIm           -       Vector of charges of the imagepoints.
! CordIm        -       Coordinates of the imagepoints.
! QImp          -       Image charge due to dipole in cavity.
! DipIm         -       Image dipole due to dipole in cavity.
! c_orbene      -       Solvent orbital energies.
! ChaNuc        -       Nuclear charges.
! qTot          -       Total charge on QM molecule.
! xyzMyQ        -       Total dipole of QM-region.
! xyzMyI        -       The induced dipole of QM-region.
! xyzMyP        -       Total dipole of the explicit solvent.
! xyzQuQ        -       Total traceless quadrupole moment of QM-region.
! CT            -       Centre of mass for QM-molecule.
!----------------------------------------------------------------------*

!----------------------------------------------------------------------*
! Common variables, unique for SCF.
!----------------------------------------------------------------------*
! INTEGER:
! --------
! iSupM         -       Pointer to the supermatrix.
! iV1           -       Pointer to MO-coefficients for QM-region.
!
! REAL:
! -----
! HHMat         -       The one-electron contribution to the Hamiltonian.
! outxyz        -       Coordinates of the MME-sites in the QM-mol.
! Cha           -       The charges in the MME expansion.
! DipMy         -       The dipoles in the MME expansion.
! Quad          -       The quadrupoles in the MME expansion.
! PotNuc        -       Nuclear repulsion.
! DenCorrD      -       Density difference between HF and MP2.
! Trace_MP2     -       Trace to MP2-HF difference density.
!----------------------------------------------------------------------*

!----------------------------------------------------------------------*
! Variables to include for the Rassi-stuff.
!----------------------------------------------------------------------*
! INTEGER:
! --------
! nState        -       Number of contracted RASSI states.
! iBigT         -       Pointer to ALL Gamma-matrices.
! nRedMO        -       Number of reduced MOs in reduced basis.
! ipAvRed       -       Pointer to optional reduced MO-basis.
!
! REAL:
! -----
! HmatState     -       The Hamiltonian matrix.
! HmatSOld      -       The stored Hamiltonian matrix.
! RasCha        -       MME-charges.
! RasDip        -       MME-dipoles.
! RasQua        -       MME-quadrupoles.
! outxyzRAS     -       The MME-centers for RASSI.
!----------------------------------------------------------------------*

!----------------------------------------------------------------------*
! iQn           -       Array that specifies the angular momentum
!                       quantum number for each basis function on
!                       solvent molecules.
! iQang         -       Like iQn, but for QM-region.
! nPrimus       -       Actually a rewriting of Icon.
! mPrimus       -       Like nPrim but for solvent.
! iCharOnBasQ   -       Charge on atom on which the ith contracted
!                       AO-basis is centered in QM-system.
! iCharOnBasC   -       Like iCharOnBasQ, but for solvent.
! iWoGehenQ     -       The (ith,jth) element tells which index the
!                       ith QM-region base (not basis-function) of
!                       the jth m_l-quantum number is to take.
!                       Needed when ordering the AO-overlaps.
! iWoGehenC     -       Like iWoGehenQ but for solvent molecule.
! nBonA_Q       -       Number of basis functions on atoms in QM.
! nBonA_C       -       Like nBonA_Q but for solvent molecule.
! CasOri        -       Array with coordinates for each basis
!                       function for solvent molecules.
! SavOri        -       Initially like CasOri, but not overwritten.
! BasOri        -       Like CasOri, but for QM-region.
! Alfa          -       Basis exponents.
! Beta          -       Like Alfa but for solvent.
! Cont          -       Contraction coefficients for QM-region.
! Dont          -       Like Cont but for solvent.
! V3            -       Original solvent MOs.
! Trans         -       Cartesian to spherical transformation.
!----------------------------------------------------------------------*

#include "maxi.fh"
integer(kind=iwp), parameter :: MxBasTri = MxBas*(MxBas+1)/2

integer(kind=iwp) :: iBigT, iCharOnBasC(MxBasC), iCharOnBasQ(MxBas), iCIInd(MxState), iCompExt(MxExtAddOn), iExtr_Atm(MxAt), &
                     iExtr_Eig, iExtra, iLuSaIn, iLuSaUt, iLuStIn, iLuStUt, iLvlShift(MxState), info_atom(MxAt), iNrExtr, iNrIn, &
                     iNrUt, Inter, iOcc1, iOrb(3), ipAvRed, iPrint, iQang(MxBas), iQn(MxBasC), iRead, iSeed, iSta, iSupM, &
                     iTcSim(64), itMax, iV1, iWoGehenC(MxBB,2*MxAngqNr-1), iWoGehenQ(MxBB,2*MxAngqNr-1), lmax, lMltSlC, &
                     mPrimus(MxBasC), nAdd, nAtom, nBA_C(MxAt), nBA_Q(MxAt), nBonA_C(3), nBonA_Q(MxAt), nCBoA_C(MxAt,MxAngqNr), &
                     nCBoA_Q(MxAt,MxAngqNr), nCent, nCha, nCIRef, nDel, nEqState, nExtAddOns, nLvlShift, nMacro, nMicro, nMlt, &
                     nPart, nPol, nPrimus(MxBas), nRedMO, NrFiles, NrStarti, NrStartu, NrStates(MxState), nSlSiteC, nState, &
                     nStateRed, nStFilT(MxParT), nTemp
real(kind=wp) :: Alfa(MxBas,MxCont), AvElcPot(MxQCen,10), BasOri(3,MxBas), Beta(MxBasC,MxCont), c_orbene(MxOrb_C), CAFieldG, &
                 CasOri(3,MxBasC), CBFieldG, CFexp, Cha(MxOT,MxQCen), ChaNuc(MxAt), CharDi(MxCen), CharDiQ(MxAt), &
                 Cont(MxBas,MxCont), CordIm(MxCen*MxPut,3), Cordst(MxCen*MxPut,3), CT(3), Cut_Elc, Cut_Ex1, Cut_Ex2, &
                 dCIRef(MxState), DelFi, DelR, DelX, DenCorrD(MxOT), Diel, DifSlExp, DipIm(MxCen*MxPut,3), DipMy(MxOT,3,MxQCen), &
                 Disp(MxPol,MxPol), dLJrep, dLvlShift(MxState), Dont(MxBasC,MxCont), Enelim, Exdt1, Exdtal, Exrep10 = Zero, &
                 Exrep2 = Zero, Exrep4 = Zero, Exrep6 = Zero, FockM(MxOT), Forcek, HHMat(MxOT), HmatSOld(MxStOT), &
                 HmatState(MxStOT) = Zero, OldGeo(MxCen*MxPut,3), outxyz(MxQCen,3), outxyzRAS(MxQCen,3), ParaTemps(MxParT), &
                 PertNElcInt(MxBasTri), Pol(MxPol), Pollim, PotNuc, Pres, QIm(MxCen*MxPut), QImp(MxCen*MxPut), Qsta(MxCha), qTot, &
                 Quad(MxOt,6,MxQCen), QuaDi(3,MxCen), QuaDiQ(3,MxAt), RasCha(MxStOT,MxQCen), RasDip(MxStOT,3,MxQCen), &
                 RasQua(MxStOT,6,MxQCen), rStart, SavOri(3,MxBasC), ScalExt(MxExtAddOn), SexRe1(MxAt,MxAt), SexRe2(MxAt,MxAt), &
                 SexRep(MxAt,MxAt), SlExpC(2,MxCen), SlExpQ(0:MxMltp,MxQCen), SlFactC(4,MxCen), SlPC(MxCen), Sqrs(MxPut*MxCen), &
                 Surf, Temp, ThrsCont, ThrsRedOcc, Trace_MP2, &
                 Trans(int(real(3*MxAngqNr**2-2*MxAngqNr-10+8*MxAngqNr**3+3*MxAngqNr**4)/real(12))), Udisp(MxAt,MxCen), &
                 V3(MxBasC,MxOrb_C), xyzMyI(3), xyzMyP(3), xyzMyQ(3), xyzQuQ(6)
logical(kind=iwp) :: AddExt, Anal, ATitle, ChargedQM, ContrStateB, DelOrAdd(12), DispDamp, EdSt, FieldDamp, lCiSelect, lExtr(12), &
                     lQuad, lSlater, MoAveRed, Mp2DensCorr, ParallelT, Qmeq, QmProd, SingPoint
character(len=100) :: JobLab
character(len=10) :: cDumpForm
character(len=8) :: ExtLabel(MxExtAddOn)
character(len=6) :: StFilIn, SaFilIn, StFilUt, SaFilUt, FieldNuc, SimEx, RassiM, EigV, QmType

public :: iTcSim, iLuStIn, iLuStUt, iLuSaIn, iLuSaUt, StFilIn, SaFilIn, StFilUt, SaFilUt, FieldNuc, SimEx, RassiM, EigV, lmax, &
          info_atom, QIm, CordIm, QImp, DipIm, c_orbene, ChaNuc, qTot, xyzMyQ, xyzMyI, xyzMyP, Sqrs, CT, xyzQuQ, iSupM, iV1, &
          HHMat, outxyz, Cha, DipMy, Quad, PotNuc, FockM, DenCorrD, Trace_MP2, nState, nStateRed, iBigT, nRedMO, ipAvRed, &
          HmatState, RasCha, RasDip, RasQua, outxyzRAS, HmatSOld, nBA_Q, nBA_C, nCBoA_Q, nCBoA_C, iQang, nPrimus, iCharOnBasQ, &
          iCharonBasC, iQn, mPrimus, iWoGehenC, iWoGehenQ, nBonA_Q, nBonA_C, Alfa, Beta, Cont, Dont, CasOri, BasOri, SavOri, V3, &
          Trans, MoAveRed, ContrStateB, lCiSelect, lExtr, lSlater, lQuad, Anal, AddExt, SingPoint, ParallelT, Mp2DensCorr, EdSt, &
          ChargedQM, DelOrAdd, ATitle, Qmeq, FieldDamp, DispDamp, QmProd, JobLab, QmType, ExtLabel, cDumpForm, AvElcPot, &
          PertNElcInt, SlFactC, SlExpC, SlPC, Cut_Elc, DifSlExp, SlExpQ, ThrsCont, dLvlShift, dCIRef, ScalExt, ParaTemps, &
          ThrsRedOcc, CFexp, Pollim, Enelim, Exdtal, Exdt1, Cut_Ex1, Cut_Ex2, CAFieldG, CBFieldG, DelX, DelFi, DelR, Forcek, &
          dLJrep, Temp, Pres, Surf, CharDi, CharDiQ, QuaDi, QuaDiQ, Udisp, Exrep2, Exrep4, Exrep6, Exrep10, Disp, Cordst, OldGeo, &
          SexRep, SexRe1, SexRe2, rStart, Diel, Qsta, Pol, nSlSiteC, lMltSlC, nMlt, iExtr_Eig, iExtr_Atm, nLvlShift, iLvlShift, &
          nCIRef, iCIInd, nExtAddOns, nTemp, nStFilT, iCompExt, NrStates, nEqState, iSeed, nMacro, nMicro, nAdd, iNrExtr, iNrIn, &
          iNrUt, NrStarti, NrStartu, nDel, NrFiles, nPol, nCha, iExtra, Inter, iPrint, iRead, iSta, itMax, iOrb, iOcc1, nPart, &
          nAtom, nCent

end module qmstat_global
