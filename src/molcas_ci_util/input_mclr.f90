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
Module Input_mclr
#include "rasdim.fh"
!
      Integer, parameter :: itociph = 64
      Integer, Parameter :: mxPer   = 255
      Integer, Parameter :: mxAtm   = MxAtom
      Integer, Parameter :: iSCF    = 1
      Integer, Parameter :: iCASSCF = 2
!
!     nactel-itoc
!
!     See RASSCF JOBIPH
!
      Integer nActEl  ! Number of active electrons
      Integer ispop   ! Spin operator
      Integer iSpin   ! Spin of state
      Integer nSym    ! Number of symmetries
      Integer State_Sym ! Symmetry of state
      Integer nConf     ! Number of configurations in state sym
      Integer nAtoms    ! Guess
      Integer mTit
      Integer lRoots    ! State from CASSCF
      Integer nRoots    ! nr of states from CASSCF
      Integer nHole1    ! nr of holes RAS1
      Integer nElec3    ! nr of el in RAS3
      Integer ntIsh,ntItri,ntIsqr,pntgrp
      Integer ntAsh,ntAtri,ntAsqr
      Integer ntBas,ntBtri,ntBsqr
      Integer nBas(mxSym)  ! Number of basis functions
      Integer nOrb(mxSym)  ! Number of basis functions
      Integer nFro(mxSym)   ! Nr of frozen orb (not allowed)
      Integer nDel(mxSym)   ! Nr of deleted orbitals (not allowed)
      Integer nOcc(mxSym)   ! nr of occ orbitals (allowed)
      Integer nExt(mxSym)   ! nr of virtual orbitals
      Integer nIsh(mxSym)   ! nr of inactive orbitals
      Integer nAsh(mxSym)   ! number of active orbitals
      Integer nRs1(mxSym)   ! number of Rs1
      Integer nRs2(mxSym)   ! nr     of rs2
      Integer nRs3(mxSym)   ! nr     of rs3
      Integer nSkip(mxSym)  ! Not active
      Integer iRoot(mxRoot) ! nonono
      Integer iPt2          ! CASPT2 orbitals
      Integer iToc(itociph)
      Integer lmode(3*mxAtm+3)  ! Which modes the response is calc for
      Integer nmode             ! Number of modes the resp is calc for
      Integer iMethod           ! 1=SCF 2=RASSCF  RELAX
      Integer niter             ! Max number of iterations
      Integer nDisp             ! Total number of perturbations
      Integer nCSF(8)           ! Number of CSFs in diff symmetries
      Integer kprint            ! print level
      Integer iBreak            ! Break criteria
!     Unit for storage of half-transformed Cho vectors
      Integer LuAChoVec(8)
!     Unit for storage of half-transformed Cho vectors
      Integer LuIChoVec(8)
      Integer LuChoInt(2)       ! Unit for storage of integrals (ii|ab)
      Integer cmsNACStates(2)   ! StateIndex NAC of MCPDFT module
      Logical isCMSNAC
!
!   nTPert
!   ======
!
!   MCKINT or rdinp
!
!   Bit 1 Singlet(0)/Triplet(1) operator
!   Bit 2 One electron contribution to perturbation (1=true)
!   Bit 3 Two electron contribution to perturbation (1=true)
!   Bit 4 Connection contribution to perturbation (1=true)
!   Bit 5 Mckinley(1) Seward(0)
!
      Integer nTPert(mxPer)
!
!     Character*8 OrbLbl(mxOrb)      ! Label on orbital ONEINT
      Character(LEN=16) Perturbation      ! The type of perturbation
      Character(LEN=LENIN) AtLbl(mxAtm) ! Label on Atom    ONEINT
      Character(LEN=72) HeaderJP(2)
      Character(LEN=72) Header1I(2)
      Character(LEN=72) TitleJp(mxTit)
      Character(LEN=80) Bline
      Character(LEN=8) TitleIN(180)
      Character(LEN=3) chirr(8)       ! Irreps
      Character(LEN=3) cmass(mxAtm)
!
!     Name of perturbation to read from ONEINT
!
      Character(LEN=8) SewLab
      Logical lCalc(3*mxAtm+3)  ! True: the perturbation
!                               !   should be calculated
      Logical ISOTOP
      Logical debug
      Logical CASINT            ! CASPT2 type integrals
      Logical page
      Logical fail              ! Calculation didnt converge
      Logical Timedep           ! Time-dependent calculation
      Logical iMCPD             ! CAS-PDFT calculation
      Logical iMSPD             ! MS-PDFT calculation
      Logical McKinley          ! Input read from McKinley
!                               !    True if MCKINT exists
      Logical TwoPert           ! Not in use
      Logical PT2               ! Read in RHS from CASPT2
      Logical SPINPOL,ELECHESS,SAVE          ! Calculate
!                                            ! spin polarization for casscf
      Logical BasCon            ! Not in use
      Logical double           ! Double isotope substitutions
      Logical multi            ! Write the response for
!                              !      diff modes to diff files
      Logical mout             ! Run only the output module -
!                              !      read resp from several files
      Logical lmass            ! Program has read user
!                              !      specified masses, umass
      Logical newCho           ! Switch to the new Cholesky algorithm
      Real*8  PotNuc,ERASSCF(mxroot),ESCF ! Energies
      Real*8  Epsilon           ! The threeshold for PCG
      Real*8  rin_ene           ! Inactive energy
      Real*8  EpsOrb            ! Not in use
      Real*8  Omega             ! Frequency of the
!                               !      time-dependent perturbation
      Real*8  Coor(3,mxAtm)
      Real*8  Weight(mxroot)
      Real*8  umass(mxAtm)      ! User specified masses
!
      Integer iEndofinput(8), imass, nfiles, nUserPT, nsRot
!
      Logical TwoStep
      Character(len=4) StepType

!
!     User-defined Temperatures and Pressure
!
      Real*8  UserP, UserT(64)
!
!     Ask Jeppe!!!
!     We all want a happy LUCIA
!
      Logical Direct
      Integer NIRREP,NSMOB,NRS0SH(1,20),NRS4SH(20,10),                  &
     &        MXR4TP, MNRS10,MXRS30
      Logical PrCI,PrOrb
      Real*8  CIthrs

       Logical RASSI

       Integer iaddressQDAT
Private :: mxRef,mxIter,mxCiIt,mxSxIt,mxTit
Private :: MaxBfn,MaxBfn_Aux,MxAO,mxAtom,mxroot,mxNemoAtom,Mxdbsc,lCache,mxact,mxina, &
           mxbas,mxOrb,mxSym,mxGAS,LENIN,LENIN1,LENIN2,LENIN3,LENIN4,LENIN5,LENIN6,LENIN8
End Module Input_mclr

