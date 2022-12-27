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
!----------------------------------------------------------------------*
! Allocate space to store the system description                       *
!----------------------------------------------------------------------*
! E1      - one-electron energy (from optimized density)               *
! E2      - two-electron energy (from optimized density)               *
! Energy  - E1 + E2                                                    *
! EKin    - kinetic energy of the system (from optimized density)      *
! E1V     - one-electron energy (variational)                          *
! E2V     - two-electron energy (variational)                          *
! EneV    - E1V + E2V                                                  *
! Elst    - E1V + E2V for each iteration, total or alpha and beta.     *
! EKinV   - kinetic energy of the system (variational)                 *
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
! ThrOcc  - threshold for occupation number to be printed              *
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
! kDisk   - array specifying position of a particular gradient         *
!           matrix in the file                                         *
! nDens   - total number of density matrices stored in the core        *
!           (note that the last matrix is just the density from        *
!           the previous iteration)                                    *
! nMem    - nDens - 1                                                  *
! nDsk    - number of density matrices stored at given instant on      *
!           disk                                                       *
! kOptim  - number of density matrices we perform optimization on      *
!           (updated for the next iteration in DIIS or INTERP)         *
! MinDMx  - maximal number of density matrices we perform minimi-      *
!           zation on (if 0 then use normal differences, maximal-      *
!           MxIter-1)                                                  *
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
! iUHF    - if = 1, perform unrestricted HF calculations               *
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
! MaxFro  - Max(nFro(i),i = 1, nSym)                                   *
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
! nOT     - Sum(nOrb(i) * (nOrb(i) + 1) / 2,i = 1, nSym)               *
! nnOc    - Sum(nOcc(i),i = 1, nSym)                                   *
! nnFr    - Sum(nFro(i),i = 1, nSym)                                   *
! nOFS    - Sum((nOrb(i)-nFro(i))**2, i = 1, nSym)                     *
! nOFT    - Sum((nOrb(i)-nFro(i))*(nOrb(i)-nFro(i)+1)/2,i=1,nSym)      *
! msymon  - Use MSYM                                                   *
!                                                                      *
! nConstr - number of constrained MOs for irrep                        *
!                                                                      *
! Atom    - name of atom *LENIN                                        *
! Type    - type of basis function *8                                  *
! Name    - Atom = Name(1:LENIN); Type = Name(LENIN+1:LENIN+8)         *
! AccCon  - acceleration convergence scheeme used in iteration         *
!                                                                      *
! ChFracMem - fraction of memory used as a buffer for keeping          *
!             Cholesky vectors previously read from disk               *
!----------------------------------------------------------------------*
Module InfSCF
use MxDM
#include "chopar.fh"

Private MxBS,MxAtms,MxIter,MxOptm,MxKeep,MxDDsk,MxTit,MxKp2U
Private MaxBfn, MaxBfn_Aux, MxAO
Private mxAtom, mxroot, mxNemoAtom, Mxdbsc, MxShll
Private mxact, mxina, mxbas, mxOrb, mxSym, mxGAS
Private LENIN, LENIN1, LENIN2, LENIN3, LENIN4, LENIN5, LENIN6, LENIN8

Integer MxConstr, nConstr(8), indxC(16,2,8), klockan
Real*8 E_nondyn


Real*8 E1,E2,EKin,PotNuc,EneV,E1V,E2V,EKinV,                      &
       EThr,DThr,FThr,DelThr,DiisTh,QudThr,ThrEne,ThrOcc,DNorm,   &
       TNorm,DMOMax,FMOMax,EDiff,Thize,RTemp,TemFac,              &
       TStop,ExFac,Tot_Charge,Tot_Nuc_Charge,Tot_El_Charge,       &
       Elst(MxIter,2),                                            &
       RotFac,RotLev,RotMax,s2uhf,HLgap,ScrFac,FlipThr,           &
       Float_End

Integer, Parameter :: nStOpt = 8
Integer :: kOptim_Max=5
Integer :: Iter_Start=1
Integer :: Iter_Ref  =1

Integer nBas(MxSym),nOrb(MxSym),nOcc(MxSym,2),                          &
        nFro(MxSym),nFrz(MxSym),nDel(MxSym),nSym,nAufb(2),              &
        isAufbauInput,nAtoms,nDens,nDsk,nMem,kOptim,                    &
        iDKeep,lPaper,MapDns(MxKeep),MapGrd(MxOptm),                    &
        iDisk(MxDDsk,2),kDisk(MxOptm),kOV(2),                           &
        nIter(0:1),nIterP,iter,jPrint,iPsLst,InVec,                     &
        kIvo,iCoCo,iUHF,jVOut,iPrOrb,iPrint,MinDMx,iDMin,               &
        MaxBas,MaxOrb,ivvloop,MaxFro,MaxOrF,MaxBxO,MaxBOO,MaxBOF,       &
        MaxOrO,nBB,nBO,nOO,nOV,mOV,nnB,nnO,nBT,nOT,nnOc,nnFr,nOFS,      &
        nOFT,nDisc,nCore,iPrForm,MaxFlip,iterprlv,                      &
        nSkip(MxSym),iAu_ab,LstVec(nStOpt),                             &
        iDummy_run, Iter2run

      Character(LEN=LENIN8), Allocatable::  Name(:)
      Character(LEN=LENIN ), Allocatable::  Atom(:)
      Character(LEN=8     ), Allocatable::  Type(:)
      Character(LEN=9     )  AccCon
      Character(LEN=4     )  Neg2_Action
      Character(LEN=80    )  KSDFT
      Character(LEN=512   )  SCF_FileOrb
      Logical isHDF5
      Integer fileorb_id
      Logical MSYMON

!----------------------------------------------------------------------*
! Allocate space for timing informations                               *
!                                                                      *
! CpuItr - CPU time per iteration                                      *
! TotCpu - total CPU time                                              *
! TotIO  - total I/O time                                              *
! TimFld - time of specified sections of the program                   *
! NamFld - names of those sections (set up in scf_init)                *
!----------------------------------------------------------------------*
      Integer,Parameter::  nFld = 16
      Real*8       CpuItr,TotCpu,TotIO,TimFld(nFld)
      Character(LEN=45) NamFld(nFld)
!----------------------------------------------------------------------*
! Allocate space to store the one-electron inntegral file header       *
!----------------------------------------------------------------------*
      Character(LEN=72) Header(2)
!----------------------------------------------------------------------*
! Allocate space to store the header of INPORB file                    *
!----------------------------------------------------------------------*
      Character(LEN=40) VTitle
      Character(LEN=80) StVec
!----------------------------------------------------------------------*
! Allocate space to store the title                                    *
!                                                                      *
! nTit  - number of title lines                                        *
! Title - title lines                                                  *
!----------------------------------------------------------------------*
      Integer nTit
      Character(LEN=72) Title(MxTit)
!----------------------------------------------------------------------*
! Allocate logical variables                                       *
!----------------------------------------------------------------------*
! DSCF       - T   = Perform direct SCF (see Subroutine OpnFls)        *
! DoCholesky - T   = Perform Cholesky SCF                              *
! DoLDF      - T   = Use Local DF fitting coefficients to construct    *
!                    Fock matrix                                       *
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
Logical DSCF,lRel,PreSch,MiniDn,WrOutD,c1Diis,Scrmbl,RFpert,Aufb,Teee,  &
        Damping,Diis,One_Grid,DoCholesky,DDnOFF,Two_Thresholds,         &
        PmTime,EmConv,WarnCfg,DoHLgap,AddFragments,WarnPocc,            &
        DoFMM,LKon,OnlyProp,DoLDF,NoExchange,WarnSlow,NoProp,           &
        FckAuf,Falcon
Logical :: RSRFO=.False., RGEK=.False.

Integer iStatPRN
End Module InfSCF
