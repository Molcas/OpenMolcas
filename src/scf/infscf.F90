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

module InfSCF

use Molcas, only: LenIn, MxSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

!----------------------------------------------------------------------*
! System description                                                   *
!----------------------------------------------------------------------*
! E1      - one-electron energy (from optimized density)               *
! E2      - two-electron energy (from optimized density)               *
! Energy  - E1 + E2                                                    *
! EKin    - kinetic energy of the system (from optimized density)      *
! E1V     - one-electron energy (variational)                          *
! E2V     - two-electron energy (variational)                          *
! EneV    - E1V + E2V                                                  *
! Elst    - E1V + E2V for each iteration, total or alpha and beta.     *
! PotNuc  - nuclear repulsion                                          *
! CPU     - CPU time per iteration                                     *
! EThr    - convergence threshold for energy                           *
! DThr    - convergence threshold for density matrix                   *
! FThr    - convergence threshold for Fock matrix                      *
! DelThr  - threshold for deleting near linear dependencies            *
! DiisTh  - threshold for density matrix to start Diis                 *
! QudThr  - threshold for Dij below which quadratic surface is         *
!           assumed                                                    *
! ThrEne  - threshold for orbital energy to be printed                 *
! DNorm   - norm of the difference of density matrices                 *
! TNorm   - corresponding norm of two-electron hamiltonian             *
! DMOMax  - the biggest off-diagonal element of density mat. in MO     *
!           basis                                                      *
! FMOMax  - the biggest off-diagonal element of Fock matrix in MO      *
!           basis                                                      *
!                                                                      *
! s2uhf   - <S^2> expectation value                                    *
!                                                                      *
! RotFac  - Scale factor for off diagonal matrix elements              *
! RotLev  - Levelshift between occupied and virtual orbitals           *
! RotMax  - Largest allowable rotation                                 *
! nSym    - number of symmetries                                       *
! nBas(i) - (i = 1, nSym), number of basis functions                   *
! nOrb(i) - (i = 1, nSym), number of orbitals                          *
! nOcc(i) - (i = 1, nSym), number of occupied orbitals                 *
! nFro(i) - (i = 1, nSym), number of frozen orbitals                   *
! nFrz(i) - (i = 1, nSym), always set to zero to compute total         *
!           density matrix                                             *
! nDel(i) - (i = 1, nSym), number of orbitals deleted by linear        *
!           dependencies                                               *
! nAufb(2)- total # occupied orbitals (if Aufbau is used...)           *
! nAtoms  - number symmetry independent of atoms in the system         *
! MapDns  - array specifying position of a particular density          *
!           matrix (positive value - core, negative - disk)            *
! iDisk   - array specifying position of a particular density          *
!           matrix in the file                                         *
! nDens   - total number of density matrices stored in the core        *
!           (note that the last matrix is just the density from        *
!           the previous iteration)                                    *
! nMem    - nDens - 1                                                  *
! kOptim  - number of density matrices we perform optimization on      *
!           (updated for the next iteration in DIIS or INTERP)         *
! iDMin   - current number of density matrices we perform minimi-      *
!           zation on                                                  *
! iDKeep  - keep-level: idkeep >2 - extended damping, extended diis    *
!                              =2 - normal   damping, extended diis    *
!                              =1 - normal   damping, normal   diis    *
!                              <1 - no       damping, no       diis    *
! nIter(0:1) - number of iterations (<= MxIter)                        *
!             nIter(0): # NDDO SCF iterations...                       *
!             nIter(1): # RHF SCF iterations:                          *
! nIterP  - either 0 or 1                                              *
! lPaper  - width of paper                                             *
! iter    - current iteration number                                   *
! iPsLst  - position of the last density                               *
! InVec   - start level                                                *
!           0 - core diagonalization,                                  *
!           1 - NDDO orbitals as intermediate step,                    *
!           2 - use INPORB,                                            *
!           3 - use DENSIN                                             *
!           4 - restart SCF calculation                                *
!           This part is under reconstruction!!!!!                     *
! LstVec  - Priority list for start of scf calculation                 *
!           This part is under reconstruction!!!!!                     *
! iCoCo   - = 1 if arbitrary occupation numbers were read              *
! kIvo    - if = 1, generate improved virtual orbitals                 *
! jVOut   - if > 1, then punch all orbitals on SCFORB                  *
! iPrOrb  - if > 0, print orbitals                                     *
! iPrForm - Format of MO output: 0-none,1-list,2-short,3-long(def)     *
! iPrint  - print level                                                *
! jPrint  - print level                                                *
! MaxBas  - Max(nBas(i),i = 1, nSym)                                   *
! MaxOrb  - Max(nOrb(i),i = 1, nSym)                                   *
! MaxOrF  - Max(nOrb(i) - nFro(i),i = 1, nSym)                         *
! MaxOrO  - Max(nOrb(i) - nOcc(i),i = 1, nSym)                         *
! MaxBxO  - Max(nBas(i) * nOrb(i),i = 1, nSym)                         *
! MaxBOF  - Max(nBas(i) * (nOrb(i) - nFro(i)),i = 1, nSym)             *
! MaxBOO  - MaxBas(n(i) * (nOrb(i) - nOcc(i)),i = 1, nSym)             *
! nBB     - Sum(nBas(i)**2,i = 1, nSym)                                *
! nOO     - Sum(nOrb(i)**2,i = 1, nSym)                                *
! nOV     - Sum((nOcc(i)-nFro(i))*(nOrb(i)-nOcc(i)),i = 1, nSym)       *
! mOV     - as nOV but without the blanks.                             *
! kOV     - as nOV but for each spin   mOV=kOV(1)+kOV(2)               *
! nBO     - Sum(nBas(i) * nOrb(i),i = 1, nSym)                         *
! nnB     - Sum(nBas(i),i = 1, nSym)                                   *
! nnO     - Sum(nOrb(i),i = 1, nSym)                                   *
! nBT     - Sum(nBas(i) * (nBas(i) + 1) / 2,i = 1, nSym)               *
! nnOc    - Sum(nOcc(i),i = 1, nSym)                                   *
! nnFr    - Sum(nFro(i),i = 1, nSym)                                   *
! nOFS    - Sum((nOrb(i)-nFro(i))**2, i = 1, nSym)                     *
! msymon  - Use MSYM                                                   *
!                                                                      *
! nConstr - number of constrained MOs for irrep                        *
!                                                                      *
! Atom    - name of atom *LenIn                                        *
! BType   - type of basis function *8                                  *
! BName   - Atom = BName(1:LenIn); BType = BName(LenIn+1:LenIn+8)      *
! AccCon  - acceleration convergence scheme used in iteration          *
!----------------------------------------------------------------------*
! For second order update info                                         *
!----------------------------------------------------------------------*
! iterso   - second order iteration number                             *
! MemRsv   - memory kept unallocated in LList management               *
! QNRTh    - threshold for QNR/C2Diis startup                          *
! DltNTh   - convergence threshold for Norm of delta                   *
! DltNrm   - actual Norm of delta after QNR/C2Diis extrapolation       *
!----------------------------------------------------------------------*
! Timing informations                                                  *
!----------------------------------------------------------------------*
! CpuItr - CPU time per iteration                                      *
! TimFld - time of specified sections of the program                   *
! NamFld - names of those sections (set up in scf_init)                *
!----------------------------------------------------------------------*
! The one-electron integral file header                                *
!----------------------------------------------------------------------*
! Header -                                                             *
!----------------------------------------------------------------------*
! The header of INPORB file                                            *
!----------------------------------------------------------------------*
! VTitle -                                                             *
! StVec  -                                                             *
!----------------------------------------------------------------------*
! The title                                                            *
!----------------------------------------------------------------------*
! nTit  - number of title lines                                        *
! Title - title lines                                                  *
!----------------------------------------------------------------------*
! Logical variables                                                    *
!----------------------------------------------------------------------*
! DSCF       - T   = Perform direct SCF (see Subroutine OpnFls)        *
! DoCholesky - T   = Perform Cholesky SCF                              *
! LKon       - T   = Perform Local-K screening in Cholesky SCF         *
! lRel       - T   = Relativistic integrals are available              *
! PreSch     - T   = Conventional prescreening                         *
! MiniDn     - T   = Use minimization of density differences           *
! DDnOFF     - T   = Use actual density as such                        *
! WrOutD     - T   = Write out last density matrix and 2el Fock matrix *
! c1Diis     - T   = perform c1-Diis acceleration                      *
! Scrmbl     - T   = add a small random nombers to core hamiltonian    *
!                    before diagonalization                            *
! RFpert     - T   = add the reaction field as a external perturbation *
! Aufb       - T   = 'Aufbau' is used...                               *
! FckAuf     - T   = Canonical orbitals occupied according to the      *
!                    Aufbau principle.                                 *
! Damping    - T   = 'Damping' is used...                              *
! Diis       - T   = 'Diis' is used...                                 *
! PmTime     - T   = Subroutine PMat_SCF is timed                      *
! EmConv     - T   = Emergency converged is signaled by someone.       *
! WarnCfg    - T   = Solution might not be lowest configuration        *
! WarnPocc   - T   = Solution with partial occ. no. found              *
! WarnSlow   - T   = Warn about slow convergence.                      *
! DoHLgap    - T   = Apply levelshift to a minimum homo-lumo gap       *
! AddFragments - F = Add the fragment's atoms to the MOLDEN file       *
!                  => to be moved to Seward due to move of GuessOrb    *
! NoExchange - Skip calculation of exchange contribution to the Fock   *
!              matrix (introduced for debugging purposes)              *
! Falcon     - T   = Fock matrix is stored in Runfile                  *
! RSRFO      - F   = Use RS-RFO instead of DIIS extrapolation          *
!----------------------------------------------------------------------*
! Globally-allocated arrays                                            *
!----------------------------------------------------------------------*
! Dens    : density matrix - vector containing some (NumDT) last       *
!           (optimized) density matrix differences - (nDT,nD,NumDT)    *
! TwoHam  : two-el. part of the Fock matrix - vector containing        *
!           corresponding 2-el. contributions - (nDT,nD,NumDT)         *
! Vxc     : Vxc     part of the Fock matrix - vector containing        *
!           corresponding 2-el. contributions - (nDT,nD,NumDT)         *
! CMO     : molecular orbitals of length nCMO                          *
! OccNo   : occupation numbers of length lthO                          *
!----------------------------------------------------------------------*

! Variables for undershot avoidance
type Loosen_Type
  real(kind=wp) :: Thrs, Thrs2, Step, Factor
end type Loosen_Type

integer(kind=iwp), parameter :: IterSO_Max = 30, MxIter = 400, MxOptm = 20, nFld = 16, nStOpt = 8

integer(kind=iwp) :: ALGO, Expand, fileorb_id, iAu_ab, iCoCo, iDisk(MxIter,2), iDKeep, iDMin, iDummy_run, indxC(16,2,8), InVec, &
                     iPrForm, iPrint, iPrOrb, iPsLst, iStatPRN, iter, Iter2run, Iter_Ref = 1, Iter_Start = 1, IterGEK = 0, &
                     iterprlv, IterSO = 0, jPrint, jVOut, kIvo, klockan, kOptim, kOptim_Max = 5, kOV(2), lPaper, LstVec(nStOpt), &
                     MapDns(MxIter), MaxBas, MaxBOF, MaxBOO, MaxBxO, MaxFlip, MaxOrb, MaxOrF, MaxOrO, MemRsv = 0, mOV, MxConstr, &
                     nAtoms, nAufb(2), nBas(MxSym), nBB, nBO, nBT, nConstr(8), nCore, nD, nDel(MxSym), nDens, nDisc, nFro(MxSym), &
                     nFrz(MxSym), nIter(0:1), nIterP, nMem, nnB, nnFr, nnO, nnOc, nOcc(MxSym,2), nOFS, nOO, nOrb(MxSym), nOV, &
                     NSCREEN, nSkip(MxSym), nSym, nTit
real(kind=wp) :: CpuItr, DE_KSDFT_c, DelThr, dFKmat, DiisTh, DltNrm = Zero, DltNTh = 0.2e-4_wp, DltNTh_old, DMOMax, &
                 dmpk = 0.045_wp, DNorm, DThr, DThr_Old, E1V, E2V, E_nondyn, EDiff, EKin, Elst(MxIter,2), Energy(MxIter) = Zero, &
                 EneV, Erest_xc, EThr, EThr_old, ExFac, FlipThr, FMOMax, FThr, FThr_Old, HLgap, PotNuc, QNRTh = 0.075_wp, QudThr, &
                 RotFac, RotLev, RotMax, RTemp, s2CNO, s2uhf, ScrFac, TemFac, Thize, ThrEne, ThrInt_old, TimFld(nFld), TNorm, &
                 Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge, TStop
logical(kind=iwp) :: AddFragments, Aufb, c1Diis, Cho_Aufb, Damping, DDnOFF, Diis, Do_Addc, Do_Tw, DoCholesky, DoFMM, DoHLgap, &
                     DSCF, EmConv, Falcon, FckAuf, isHDF5, LKon, lRel, MiniDn, MSYMON, NoExchange, NoFerm, NoProp, One_Grid, &
                     OnlyProp, PmTime, PreSch, REORD, RFpert, RGEK, RSRFO, Scrmbl, Teee, Two_Thresholds, WarnCfg, WarnPocc, &
                     WarnSlow, WrOutD
character(len=512) :: SCF_FileOrb
character(len=80) :: ADDC_KSDFT, KSDFT, StVec
character(len=72) :: Header(2), Title(1)
character(len=45) :: NamFld(nFld)
character(len=40) :: VTitle
character(len=9) :: AccCon
character(len=4) :: Neg2_Action
type(Loosen_Type) :: Loosen
integer(kind=iwp), allocatable :: OrbType(:,:)
real(kind=wp), allocatable :: CMO(:,:), CMO_ref(:,:), Darwin(:), EDFT(:), EOrb(:,:), FockAO(:,:), HDiag(:), KntE(:), MssVlc(:), &
                              OccNo(:,:), OneHam(:), Ovrlp(:), TrDD(:,:,:), TrDh(:,:,:), TrDP(:,:,:), TrM(:,:)
real(kind=wp), allocatable :: OccSet_e(:,:), OccSet_m(:,:)
real(kind=wp), allocatable, target :: Dens(:,:,:), FockMO(:,:), TwoHam(:,:,:), Vxc(:,:,:)
character(len=LenIn+8), allocatable :: BName(:)
character(len=LenIn), allocatable :: Atom(:)
character(len=8), allocatable :: BType(:)

public :: AccCon, ADDC_KSDFT, AddFragments, ALGO, Atom, Aufb, BName, BType, c1Diis, Cho_Aufb, CMO, CMO_ref, CpuItr, Damping, &
          Darwin, DDnOFF, DE_KSDFT_c, DelThr, Dens, dFKmat, Diis, DiisTh, DltNrm, DltNth, DltNTh_old, DMOMax, dmpk, DNorm, &
          Do_Addc, Do_Tw, DoCholesky, DoFMM, DoHLgap, DSCF, DThr, DThr_Old, E1V, E2V, E_nondyn, EDFT, EDiff, EKin, Elst, EmConv, &
          Energy, EneV, EOrb, Erest_xc, EThr, EThr_old, ExFac, Expand, Falcon, FckAuf, fileorb_id, FlipThr, FMOMax, FockAO, &
          FockMO, FThr, FThr_Old, HDiag, Header, HLgap, iAu_ab, iCoCo, iDisk, iDKeep, iDMin, iDummy_run, indxC, InVec, iPrForm, &
          iPrint, iPrOrb, iPsLst, isHDF5, iStatPRN, iter, Iter2run, Iter_Ref, Iter_Start, IterGEK, iterprlv, IterSO, IterSO_Max, &
          jPrint, jVOut, kIvo, klockan, KntE, kOptim, kOptim_Max, kOV, KSDFT, LKon, Loosen, lPaper, lRel, LstVec, MapDns, MaxBas, &
          MaxBOF, MaxBOO, MaxBxO, MaxFlip, MaxOrb, MaxOrF, MaxOrO, MemRsv, MiniDn, mOV, MssVlc, MSYMON, MxConstr, MxIter, MxOptm, &
          NamFld, nAtoms, nAufb, nBas, nBB, nBO, nBT, nConstr, nCore, nD, nDel, nDens, nDisc, Neg2_Action, nFld, nFro, nFrz, &
          nIter, nIterP, nMem, nnB, nnFr, nnO, nnOc, nOcc, NoExchange, NoFerm, nOFS, nOO, NoProp, nOrb, nOV, NSCREEN, nSkip, &
          nStOpt, nSym, nTit, OccNo, OccSet_e, OccSet_m, One_Grid, OneHam, OnlyProp, OrbType, Ovrlp, PmTime, PotNuc, PreSch, &
          QNRTh, QudThr, REORD, RFpert, RGEK, RotFac, RotLev, RotMax, RSRFO, RTemp, s2CNO, s2uhf, SCF_FileOrb, ScrFac, Scrmbl, &
          StVec, Teee, TemFac, Thize, ThrEne, ThrInt_old, TimFld, Title, TNorm, Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge, TrDD, &
          TrDh, TrDP, TrM, TStop, Two_Thresholds, TwoHam, VTitle, Vxc, WarnCfg, WarnPocc, WarnSlow, WrOutD

end module InfSCF
