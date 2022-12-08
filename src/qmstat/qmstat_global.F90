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

use Definitions, only: wp, iwp

implicit none
private

!----------------------------------------------------------------------*
! Define some bounds for qmstat.
!----------------------------------------------------------------------*
! MxAngqNr      -        Maximal angular quantum number for basis
!                        functions in QM-region. 1 is s, 2 is p etc.
! MxMltp        -        Highest order multipole in MME
! MxPut         -        Maximal number of molecules to put in system
! MxSymQ        -        Maximal number of symmetries
!----------------------------------------------------------------------*

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
! SupM          -       The supermatrix.
! V1            -       MO-coefficients for QM-region.
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
! BigT          -       ALL Gamma-matrices.
! nRedMO        -       Number of reduced MOs in reduced basis.
! AvRed         -       Optional reduced MO-basis.
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

integer(kind=iwp), parameter :: MxAngqNr = 7, MxMltp = 3, MxPut = 220, MxSymQ = 1

integer(kind=iwp) :: iExtr_Eig, iExtra, iLuSaIn, iLuSaUt, iLuStIn, iLuStUt, iNrIn, iNrUt, Inter, iOcc1, iOrb(3), iPrint, iRead, &
                     iSeed, iSta, iTcSim(64), itMax, lmax, lMltSlC, nAdd, nAtom, nCent, nCha, nCIRef, nDel, nEqState, nExtAddOns, &
                     nLvlShift, nMacro, nMicro, nMlt, nPart, nPol, nRedMO, NrFiles, NrStarti, NrStartu, nSlSiteC, nState, nTemp
real(kind=wp) :: CAFieldG, CBFieldG, CFexp, CharDi(2), CT(3), Cut_Elc, Cut_Ex1, Cut_Ex2, DelFi, DelR, DelX, Diel, DifSlExp, &
                 dLJrep, Enelim, Exdt1, Exdtal, Exrep10, Exrep2, Exrep4, Exrep6, Forcek, Pollim, PotNuc, Pres, qTot, QuaDi(3,2), &
                 rStart, Surf, Temp, ThrsCont, ThrsRedOcc, Trace_MP2, xyzMyI(3), xyzMyP(3), xyzMyQ(3), xyzQuQ(6)
logical(kind=iwp) :: AddExt, Anal, ATitle, ChargedQM, ContrStateB, DelOrAdd(12), DispDamp, EdSt, FieldDamp, lCiSelect, lExtr(12), &
                     lQuad, lSlater, MoAveRed, Mp2DensCorr, ParallelT, Qmeq, QmProd, SingPoint
character(len=100) :: JobLab
character(len=10) :: cDumpForm
character(len=6) :: StFilIn, SaFilIn, StFilUt, SaFilUt, FieldNuc, SimEx, RassiM, EigV, QmType
integer(kind=iwp), allocatable :: iCIInd(:), iCompExt(:), iExtr_Atm(:), iLvlShift(:), info_atom(:), iQang(:), iQn(:), &
                                  iWoGehenC(:,:), iWoGehenQ(:,:), mPrimus(:), nPrimus(:), nBA_C(:), nBA_Q(:), nBonA_C(:), &
                                  nBonA_Q(:), nCBoA_C(:,:), nCBoA_Q(:,:), nCnC_C(:), NrStates(:), nStFilT(:)
real(kind=wp), allocatable :: Alfa(:,:), AvElcPot(:,:), AvRed(:,:), BasOri(:,:), Beta(:,:), BigT(:,:), c_orbene(:), CasOri(:,:), &
                              Cha(:,:), ChaNuc(:), CharDiQ(:), Cont(:,:), CordIm(:,:), Cordst(:,:), dCIRef(:), DenCorrD(:), &
                              DipIm(:,:), DipMy(:,:,:), Disp(:,:), dLvlShift(:), Dont(:,:), FockM(:), HHmat(:), HmatSOld(:), &
                              HmatState(:), OldGeo(:,:), outxyz(:,:), outxyzRAS(:,:), Paratemps(:), PertNElcInt(:), Pol(:), &
                              QIm(:), QImp(:), Qsta(:), Quad(:,:,:), QuaDiQ(:,:), RasCha(:,:), RasDip(:,:,:), RasQua(:,:,:), &
                              SavOri(:,:), ScalExt(:), Sexre1(:,:), Sexre2(:,:), Sexrep(:,:), SlExpC(:,:), SlExpQ(:,:), &
                              SlFactC(:,:), SlPC(:), Sqrs(:), SupM(:,:), Trans(:), Udisp(:,:), V1(:,:), V3(:,:)
character(len=8), allocatable :: ExtLabel(:)

public :: AddExt, Alfa, Anal, ATitle, AvElcPot, AvRed, BasOri, Beta, BigT, c_orbene, CAFieldG, CasOri, CBFieldG, cDumpForm, CFexp, &
          Cha, ChaNuc, CharDi, CharDiQ, ChargedQM, Cont, ContrStateB, CordIm, Cordst, CT, Cut_Elc, Cut_Ex1, Cut_Ex2, dCIRef, &
          DelFi, DelOrAdd, DelR, DelX, DenCorrD, Diel, DifSlExp, DipIm, DipMy, Disp, DispDamp, dLJrep, dLvlShift, Dont, EdSt, &
          EigV, Enelim, Exdt1, Exdtal, Exrep10, Exrep2, Exrep4, Exrep6, ExtLabel, FieldDamp, FieldNuc, FockM, Forcek, HHMat, &
          HmatSOld, HmatState, iCIInd, iCompExt, iExtr_Atm, iExtr_Eig, iExtra, iLuSaIn, iLuSaUt, iLuStIn, iLuStUt, iLvlShift, &
          info_atom, iNrIn, iNrUt, Inter, iOcc1, iOrb, iPrint, iQang, iQn, iRead, iSeed, iSta, iTcSim, itMax, iWoGehenC, &
          iWoGehenQ, JobLab, lCiSelect, lExtr, lmax, lMltSlC, lQuad, lSlater, MoAveRed, Mp2DensCorr, mPrimus, MxAngqNr, MxMltp, &
          MxPut, MxSymQ, nAdd, nAtom, nBA_C, nBA_Q, nBonA_C, nBonA_Q, nCBoA_C, nCBoA_Q, nCent, nCha, nCIRef, nCnC_C, nDel, &
          nEqState, nExtAddOns, nLvlShift, nMacro, nMicro, nMlt, nPart, nPol, nPrimus, nRedMO, NrFiles, NrStarti, NrStartu, &
          NrStates, nSlSiteC, nState, nStFilT, nTemp, OldGeo, outxyz, outxyzRAS, ParallelT, ParaTemps, PertNElcInt, Pol, Pollim, &
          PotNuc, Pres, QIm, QImp, Qmeq, QmProd, Qmstat_end, QmType, Qsta, qTot, Quad, QuaDi, QuaDiQ, RasCha, RasDip, RasQua, &
          RassiM, rStart, SaFilIn, SaFilUt, SavOri, ScalExt, SexRe1, SexRe2, SexRep, SimEx, SingPoint, SlExpC, SlExpQ, SlFactC, &
          SlPC, Sqrs, StFilIn, StFilUt, SupM, Surf, Temp, ThrsCont, ThrsRedOcc, Trace_MP2, Trans, Udisp, V1, V3, xyzMyI, xyzMyP, &
          xyzMyQ, xyzQuQ

contains

subroutine Qmstat_end()

  use stdalloc, only: mma_deallocate

  if (lCiSelect) then
    call mma_deallocate(iCIInd)
    call mma_deallocate(dCIRef)
  end if

  if (nLvlShift > 0) then
    call mma_deallocate(iLvlShift)
    call mma_deallocate(dLvlShift)
  end if

  if (AddExt) then
    call mma_deallocate(ScalExt)
    call mma_deallocate(ExtLabel)
    call mma_deallocate(iCompExt)
  end if

  if (ParallelT) then
    call mma_deallocate(nStFilT)
  end if

  if (DispDamp) then
    call mma_deallocate(CharDiQ)
    call mma_deallocate(QuaDiQ)
  end if

  if (allocated(Alfa)) call mma_deallocate(Alfa)
  if (allocated(AvElcPot)) call mma_deallocate(AvElcPot)
  if (allocated(Beta)) call mma_deallocate(Beta)
  if (allocated(c_orbene)) call mma_deallocate(c_orbene)
  if (allocated(Cha)) call mma_deallocate(Cha)
  if (allocated(Cont)) call mma_deallocate(Cont)
  if (allocated(CordIm)) call mma_deallocate(CordIm)
  if (allocated(DenCorrD)) call mma_deallocate(DenCorrD)
  if (allocated(DipMy)) call mma_deallocate(DipMy)
  if (allocated(Dont)) call mma_deallocate(Dont)
  if (allocated(HHmat)) call mma_deallocate(HHmat)
  if (allocated(HmatSOld)) call mma_deallocate(HmatSOld)
  if (allocated(HmatState)) call mma_deallocate(HmatState)
  if (allocated(NrStates)) call mma_deallocate(NrStates)
  if (allocated(outxyz)) call mma_deallocate(outxyz)
  if (allocated(outxyzRAS)) call mma_deallocate(outxyzRAS)
  if (allocated(Quad)) call mma_deallocate(Quad)
  if (allocated(RasCha)) call mma_deallocate(RasCha)
  if (allocated(RasDip)) call mma_deallocate(RasDip)
  if (allocated(RasQua)) call mma_deallocate(RasQua)
  if (allocated(SlExpQ)) call mma_deallocate(SlExpQ)
  if (allocated(Udisp)) call mma_deallocate(Udisp)
  if (allocated(V3)) call mma_deallocate(V3)

  call mma_deallocate(Cordst)
  call mma_deallocate(Disp)
  call mma_deallocate(iExtr_Atm)
  call mma_deallocate(Pol)
  call mma_deallocate(Qsta)
  call mma_deallocate(Sexre1)
  call mma_deallocate(Sexre2)
  call mma_deallocate(Sexrep)
  call mma_deallocate(SlExpC)
  call mma_deallocate(SlFactC)
  call mma_deallocate(SlPC)

end subroutine Qmstat_end

end module qmstat_global
