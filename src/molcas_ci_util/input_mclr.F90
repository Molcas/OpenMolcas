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

module Input_MCLR

! nActEl                : Number of active electrons
! iSpin                 : Spin of state
! nSym                  : Number of symmetries
! State_Sym             : Symmetry of state
! nConf                 : Number of configurations in state sym
! nAtoms                : Guess
! lRoots                : State from CASSCF
! nRoots                : nr of states from CASSCF
! nHole1                : nr of holes RAS1
! nElec3                : nr of el in RAS3
! nBas                  : Number of basis functions
! nOrb                  : Number of basis functions
! nFro                  : Nr of frozen orb (not allowed)
! nDel                  : Nr of deleted orbitals (not allowed)
! nIsh                  : nr of inactive orbitals
! nAsh                  : number of active orbitals
! nRs1                  : number of Rs1
! nRs2                  : nr     of rs2
! nRs3                  : nr     of rs3
! nSkip                 : Not active
! iRoot                 : nonono
! iPt2                  : CASPT2 orbitals
! iMethod               : 1=SCF 2=RASSCF  RELAX
! niter                 : Max number of iterations
! nDisp                 : Total number of perturbations
! nCSF                  : Number of CSFs in diff symmetries
! kprint                : print level
! iBreak                : Break criteria
! LuAChoVec             : Unit for storage of half-transformed Cho vectors
! LuIChoVec             : Unit for storage of half-transformed Cho vectors
! LuChoInt              : Unit for storage of integrals (ii|ab)
! Perturbation          : The type of perturbation
! AtLbl                 : Label on Atom    ONEINT
! chirr                 : Irreps
! lCalc                 : True: the perturbation should be calculated
! CASINT                : CASPT2 type integrals
! fail                  : Calculation didnt converge
! Timedep               : Time-dependent calculation
! iMCPD                 : CAS-PDFT calculation
! iMSPD                 : MS-PDFT calculation
! McKinley              : Input read from McKinley True if MCKINT exists
! PT2                   : Read in RHS from CASPT2
! SPINPOL               : Calculate spin polarization for casscf
! double                : Double isotope substitutions
! newCho                : Switch to the new Cholesky algorithm
! PotNuc, ERASSCF, ESCF : Energies
! Eps                   : The threshold for PCG
! rin_ene               : Inactive energy
! Omega                 : Frequency of the time-dependent perturbation
! UserP, UserT          : User-defined Temperatures and Pressure
!
! nTPert
! ======
!
! MCKINT or rdinp
!
! Bit 1 Singlet(0)/Triplet(1) operator
! Bit 2 One electron contribution to perturbation (1=true)
! Bit 3 Two electron contribution to perturbation (1=true)
! Bit 4 Connection contribution to perturbation (1=true)
! Bit 5 McKinley(1) Seward(0)

use Molcas, only: LenIn, MxAtom, MxRoot, MxSym
use RASDim, only: MxTit
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: itociph = 64, mxPer = 255

integer(kind=iwp) :: iaddressQDAT, iBreak, iMethod, iPt2, iRoot(mxRoot), iSpin, ispop, iToc(itociph), kprint, lRoots, &
                     LuAChoVec(8), LuChoInt(2), LuIChoVec(8), mTit, nActEl, nAsh(mxSym), nAtoms, nBas(mxSym), nConf, nCSF(8), &
                     nDel(mxSym), nDisp, nElec3, nFro(mxSym), nHole1, NIRREP, nIsh(mxSym), niter, nOrb(mxSym), nRoots, &
                     nRs1(mxSym), nRs2(mxSym), nRs3(mxSym), nSkip(mxSym), nsRot, nSym, ntAsh, ntAsqr, ntAtri, ntBas, ntBsqr, &
                     ntBtri, ntIsh, ntIsqr, ntItri, nTPert(mxPer), nUserPT, State_Sym
real(kind=wp) :: Coor(3,MxAtom), Eps, ERASSCF(mxroot), ESCF, Omega, PotNuc, rin_ene, UserP, UserT(64), Weight(mxroot)
logical(kind=iwp) :: CASINT, debug, double, fail, iMCPD, iMSPD, lCalc(3*MxAtom+3), McKinley, newCho, page, PT2, RASSI, lSave, &
                     SPINPOL, Timedep, TwoStep
character(len=LenIn) :: AtLbl(MxAtom)
character(len=72) :: Header1I(2), HeaderJP(2), TitleJp(mxTit)
character(len=16) :: Perturbation
character(len=8) :: TitleIN(180)
character(len=4) :: StepType
character(len=3) :: chirr(8)

public :: AtLbl, CASINT, chirr, Coor, debug, double, Eps, ERASSCF, ESCF, fail, Header1I, HeaderJP, iaddressQDAT, iBreak, iMCPD, &
          iMethod, iMSPD, iPt2, iRoot, iSpin, ispop, iToc, itociph, kprint, lCalc, lRoots, LuAChoVec, LuChoInt, LuIChoVec, &
          McKinley, mTit, nActEl, nAsh, nAtoms, nBas, nConf, nCSF, nDel, nDisp, nElec3, newCho, nFro, nHole1, NIRREP, nIsh, niter, &
          nOrb, nRoots, nRs1, nRs2, nRs3, nSkip, nsRot, nSym, ntAsh, ntAsqr, ntAtri, ntBas, ntBsqr, ntBtri, ntIsh, ntIsqr, ntItri, &
          nTPert, nUserPT, Omega, page, Perturbation, PotNuc, PT2, RASSI, rin_ene, lSave, SPINPOL, State_Sym, StepType, Timedep, &
          TitleIN, TitleJp, TwoStep, UserP, UserT, Weight

end module Input_MCLR
