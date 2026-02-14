!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Naoki Nakatani                                         *
!***********************************************************************
!  DMRGCtl
!
!> @brief
!>   DMRG Control
!> @author N. Nakatani
!>
!> @details
!> Many functionality is same as ::CICtl (copy most part of ::CICtl to here)
!> To use DMRG as CAS-CI solver.
!> \p IRst = ``0`` is used for the first and the last (DMRG-) CI calculations
!> to fully optimize the DMRG wavefunction.
!> \p IRst = ``1`` enables to restart DMRG calculation from last iteration
!>
!> @param[in]     CMO    MO coefficients
!> @param[out]    D      Average 1-dens matrix
!> @param[out]    DS     Average spin 1-dens matrix
!> @param[out]    P      Average symm. 2-dens matrix
!> @param[out]    PA     Average antisymm. 2-dens matrix
!> @param[out]    FI     Fock matrix from inactive density
!> @param[in,out] D1I    Inactive 1-dens matrix
!> @param[in,out] D1A    Active 1-dens matrix
!> @param[in]     TUVX   Active 2-el integrals
!> @param[in]     IFINAL Calculation status switch
!> @param[in]     IRst   DMRG restart status switch
!***********************************************************************

#include "compiler_features.h"

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
subroutine DMRGCtl(CMO,D,DS,P,PA,FI,D1I,D1A,TUVX,IFINAL,IRst)

use wadr, only: FMO
use rctfld_module, only: lRF
use casvb_global, only: ifvb
use timers, only: TimeDens
use lucia_data, only: DStmp, Dtmp, PAtmp, Pscr, Ptmp
use gas_data, only: iDoGAS
use rasscf_global, only: DFTFOCK, ExFac, iAdr15, iPCMRoot, iRoot, ITER, KSDFT, lRoots, n_Det, NAC, NACPAR, NACPR2, nFint, nRoots, &
                         S, Weight
use PrintLevel, only: DEBUG, INSANE
use output_ras, only: IPRLOC
use general_data, only: ISPIN, jobiph, nactel, nash, ntot2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMO(*), D(*), DS(*), P(*), PA(*), FI(*), D1I(*), D1A(*), TUVX(*)
integer(kind=iwp) :: iFinal, IRst
integer(kind=iwp) :: i, iPrLev, jDisk, jRoot, kRoot, NACT4, nTmpPUVX
real(kind=wp) :: dum1, dum2, dum3, rdum(1), Scal, Time(2)
logical(kind=iwp) :: Do_ESPF
real(kind=wp), allocatable :: P2MO(:), RCT(:), RCT_F(:), RCT_FS(:), RCT_S(:), TmpD1S(:), TmpDS(:), TmpPUVX(:), TmpTUVX(:)

IPRLEV = IPRLOC(3)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering DMRGCTL'

! set up flag 'IFCAS' for GAS option, which is set up in gugatcl originally.
! IFCAS = 0: This is a CAS calculation
! IFCAS = 1: This is a RAS calculation - This might cause an error in DMRG-CASSCF.

if (iDoGas) call setsxci()
if (IPRLEV > DEBUG) then
  write(u6,*)
  write(u6,*) ' Enter DMRG section'
  write(u6,*) ' =================='
  write(u6,*)
  write(u6,*) ' iteration count =',ITER
end if

! SOME DIRTY SETUPS

S = Half*real(ISPIN-1,kind=wp)

! COMPUTE ONE ELECTRON INTEGRALS IN MO BASIS
! AND ADD CORE INTERACTION

! FMO FOCK MATRIX IN MO-BASIS
! LW2: 1-PARTICLE DENSITY MATRIX ALSO USED IN MO/AO TRANSFORMATION

call mma_allocate(FMO,NACPAR,Label='FMO')
call DecideOnESPF(Do_ESPF)
if (lRF .or. (KSDFT /= 'SCF') .or. Do_ESPF) then

  ! In case of a reaction field in combination with an average CAS
  ! select the potential of the appropriate state.

  jDisk = IADR15(3)
  do i=1,IPCMROOT-1
    call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
    call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
    call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
    call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
  end do

  call mma_allocate(RCT_F,NTOT2,Label='RCT_F')
  call mma_allocate(RCT_FS,NTOT2,Label='RCT_FS')
  if (IFINAL == 0) then

    ! Use normal MOs

    call mma_allocate(RCT,NACPAR,Label='RCT')
    call mma_allocate(P2MO,NACPR2,Label='P2MO')

    ! Get the total density in MOs

    call DDafile(JOBIPH,2,RCT,NACPAR,jDisk)
    call Put_dArray('D1mo',RCT,NACPAR)  ! Put it on the RUNFILE
    if (NASH(1) /= NAC) call DBLOCK(RCT)
    ! Transform to AOs
    call Get_D1A_RASSCF(CMO,RCT,RCT_F)

    ! Get the spin density in MOs

    if (NACTEL == 0) then
      call DCOPY_(NTOT2,[Zero],0,RCT_FS,1)
    else
      call mma_allocate(RCT_S,NACPAR,Label='RCT_S')
      call DDafile(JOBIPH,2,RCT_S,NACPAR,jDisk)
      if (NASH(1) /= NAC) call DBLOCK(RCT_S)
      ! Transform to AOs
      call Get_D1A_RASSCF(CMO,RCT_S,RCT_FS)
      call mma_deallocate(RCT_S)
    end if

    ! Get the 2-particle density in MO

    call DDafile(JOBIPH,2,P2MO,NACPR2,jDisk)
    call Put_dArray('P2mo',P2MO,NACPR2) ! Put it on the RUNFILE

    call SGFCIN(CMO,FMO,FI,D1I,RCT_F,RCT_FS)

    call mma_deallocate(P2MO)
    call mma_deallocate(RCT)

  else

    ! Here the pseudo-natural orbitals are in CMO and we need to
    ! get the D1A of the selected state in this basis.

    ! Compute the density of the particular state

    call mma_allocate(Dtmp,NAC**2,Label='Dtmp')
    call mma_allocate(DStmp,NAC**2,Label='DStmp')
    call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
    if (NAC >= 1) then
      if (NACTEL == 0) then
        Dtmp(:) = Zero
        DStmp(:) = Zero
        Ptmp(:) = Zero
      else
        ! load back 1- and 2-RDMs from previous DMRG run
        NACT4 = NAC**4
        call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
        call mma_allocate(Pscr,NACT4,Label='Pscr')
#       ifdef _ENABLE_BLOCK_DMRG_
        call block_densi_rasscf(IPCMRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#       elif _ENABLE_CHEMPS2_DMRG_
        call chemps2_densi_rasscf(IPCMRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#       elif _ENABLE_DICE_SHCI_
        call dice_densi_rasscf(IPCMRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#       endif

        ! NN.14 NOTE: IFCAS must be 0 for DMRG-CASSCF
        !if (IFCAS > 2) call CISX(IDXSX,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
        call mma_deallocate(Pscr)
        call mma_deallocate(PAtmp)
      end if

    else
      Dtmp(:) = Zero
      DStmp(:) = Zero
      Ptmp(:) = Zero
    end if
    ! Modify the symmetric 2-particle density if only partial
    ! "exact exchange" is included.
    !n_Det = 2
    !n_unpaired_elec = iSpin-1
    !n_paired_elec = nActEl-n_unpaired_elec
    !if (n_unpaired_elec+n_paired_elec/2 == nac) n_Det = 1
    if (ExFac /= One) call Mod_P2(Ptmp,NACPR2,Dtmp,NACPAR,DStmp,ExFac,n_Det)

    call Put_dArray('P2mo',Ptmp,NACPR2) ! Put it on the RUNFILE

    call mma_deallocate(Ptmp)

    call Put_dArray('D1mo',Dtmp,NACPAR) ! Put it on the RUNFILE
    if (NASH(1) /= NAC) call DBLOCK(Dtmp)
    call Get_D1A_RASSCF(CMO,Dtmp,RCT_F)

    if (NASH(1) /= NAC) call DBLOCK(DStmp)
    call Get_D1A_RASSCF(CMO,DStmp,RCT_FS)

    call mma_deallocate(Dtmp)
    call mma_deallocate(DStmp)

    call SGFCIN(CMO,FMO,FI,D1I,RCT_F,RCT_FS)

  end if
  call mma_deallocate(RCT_FS)
  call mma_deallocate(RCT_F)

else
  ! Normal case

  call mma_allocate(TmpDS,NACPAR,Label='TmpDS')
  call mma_allocate(TmpD1S,NTOT2,Label='TmpD1S')
  call dcopy_(NACPAR,DS,1,TmpDS,1)
  if (NASH(1) /= NAC) call DBLOCK(TmpDS)
  call Get_D1A_RASSCF(CMO,TmpDS,TmpD1S)
  call mma_deallocate(TmpDS)

  call SGFCIN(CMO,FMO,FI,D1I,D1A,TmpD1S)
  call mma_deallocate(TmpD1S)

end if

if (IfVB == 2) goto 9000

! SOLVE DMRG WAVEFUNCTION

if (IfVB == 1) then
  ! NN.14 FIXME: I'm not sure whether this option should work?
  call cvbmn_rvb(max(ifinal,1))
else
  if ((KSDFT(1:3) /= 'SCF') .and. (DFTFOCK(1:4) == 'DIFF') .and. (nac /= 0)) then
    nTmpPUVX = nFint
    call mma_allocate(TmpPUVX,nTmpPUVX,Label='TmpPUVX')
    call mma_allocate(TmpTUVX,NACPR2,Label='TmpTUVX')
    TmpTUVX(:) = Zero
    call Get_dArray('DFT_TwoEl',TmpPUVX,nTmpPUVX)
    call Get_TUVX(TmpPUVX,TmpTUVX)
    call DaXpY_(NACPR2,One,TUVX,1,TmpTUVX,1)
#   ifdef _ENABLE_BLOCK_DMRG_
    call BlockCtl(FMO,TmpTUVX,IFINAL,IRst)
#   elif _ENABLE_CHEMPS2_DMRG_
    call Chemps2Ctl(FMO,TmpTUVX,IFINAL,IRst)
#   elif _ENABLE_DICE_SHCI_
    call DiceCtl(FMO,TmpTUVX,IFINAL,IRst)
#   endif

    call mma_deallocate(TmpTUVX)
    call mma_deallocate(TmpPUVX)
  else
#   ifdef _ENABLE_BLOCK_DMRG_
    call BlockCtl(FMO,TUVX,IFINAL,IRst)
#   elif _ENABLE_CHEMPS2_DMRG_
    call Chemps2Ctl(FMO,TUVX,IFINAL,IRst)
#   elif _ENABLE_DICE_SHCI_
    call DiceCtl(FMO,TUVX,IFINAL,IRst)
#   endif
  end if
end if

! CALCULATE DENSITY MATRICES
! SAVE DENSITY MATRICES ON FILE
! COMPUTE AVERAGE DENSITY MATRICES

! Dtmp: ONE-BODY DENSITY
! DStmp: ONE-BODY SPIN DENSITY
! Ptmp: SYMMETRIC TWO-BODY DENSITY
! PAtmp: ANTISYMMETRIC TWO-BODY DENSITY

call Timing(Time(1),dum1,dum2,dum3)
call dCopy_(NACPAR,[Zero],0,D,1)
call dCopy_(NACPAR,[Zero],0,DS,1)
call dCopy_(NACPR2,[Zero],0,P,1)
call dCopy_(NACPR2,[Zero],0,PA,1)
call mma_allocate(Dtmp,NAC**2,Label='Dtmp')
call mma_allocate(DStmp,NAC**2,Label='DStmp')
call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
jDisk = IADR15(3)

do jRoot=1,lRoots
  ! load density matrices from DMRG run
  if (NAC >= 1) then
    NACT4 = NAC**4
    call mma_allocate(Pscr,NACT4,Label='Pscr')
#   ifdef _ENABLE_BLOCK_DMRG_
    call block_densi_rasscf(jRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#   elif _ENABLE_CHEMPS2_DMRG_
    call chemps2_densi_rasscf(jRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#   elif _ENABLE_DICE_SHCI_
    call dice_densi_rasscf(jRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#   endif
    call mma_deallocate(Pscr)
  end if
  ! Modify the symmetric 2-particle density if only partial
  ! "exact exchange" is included.
  !n_Det = 2
  !n_unpaired_elec = iSpin-1
  !n_paired_elec = nActEl-n_unpaired_elec
  !if (n_unpaired_elec+n_paired_elec/2 == nac) n_Det = 1
  !  write(u6,*) ' iSpin=',iSpin
  !  write(u6,*) ' n_unpaired_elec',n_unpaired_elec
  !  write(u6,*) ' n_paired_elec',n_paired_elec
  !  write(u6,*) ' n_unpaired_elec+n_paired_elec/2',n_unpaired_elec+n_paired_elec/2
  !  write(u6,*) ' n_Det=',n_Det
  !end if

  if (ExFac /= One) call Mod_P2(Ptmp,NACPR2,Dtmp,NACPAR,DStmp,ExFac,n_Det)

  ! update average density matrices
  Scal = Zero
  do kRoot=1,nRoots
    if (iRoot(kRoot) == jRoot) then
      Scal = Weight(kRoot)
      exit
    end if
  end do
  call daxpy_(NACPAR,Scal,Dtmp,1,D,1)
  call daxpy_(NACPAR,Scal,DStmp,1,DS,1)
  call daxpy_(NACPR2,Scal,Ptmp,1,P,1)
  call daxpy_(NACPR2,Scal,PAtmp,1,PA,1)
  ! save density matrices on disk
  call DDafile(JOBIPH,1,Dtmp,NACPAR,jDisk)
  call DDafile(JOBIPH,1,DStmp,NACPAR,jDisk)
  call DDafile(JOBIPH,1,Ptmp,NACPR2,jDisk)
  call DDafile(JOBIPH,1,PAtmp,NACPR2,jDisk)
end do

call mma_deallocate(PAtmp)
call mma_deallocate(Ptmp)
call mma_deallocate(DStmp)
call mma_deallocate(Dtmp)

! PREPARE DENSITY MATRICES AS USED BY THE SUPER CI SECTION

! print matrices
if (IPRLEV >= INSANE) then
  call TRIPRT('Averaged one-body density matrix, D',' ',D,NAC)
  call TRIPRT('Averaged one-body spin density matrix, DS',' ',DS,NAC)
  call TRIPRT('Averaged two-body density matrix, P',' ',P,NACPAR)
  call TRIPRT('Averaged antisymmetric two-body density matrix,PA',' ',PA,NACPAR)
end if
if (NASH(1) /= NAC) call DBLOCK(D)
call Timing(Time(2),dum1,dum2,dum3)
TimeDens = TimeDens+Time(2)-Time(1)

call mma_deallocate(FMO)

9000 continue

! For RF calculations make sure that the we are following the
! correct root.

! In the current implementation the overlap between the CI vectors
! of different macro iterations is used. This criterion stricktly
! only hold if the orbitals are not changed in between the
! interations, to make sure that this approximately holds the
! comparision is only to be considered to be valid if the rotmax
! parameter is below an empirical threshold. In the future the
! procedure for indentifying root flipping has to be made more
! robust.
!
!SVC: if CISElect is used, roots should be traced regardless of orbital
!     rotation, so in that case ignore the automatic tracing and follow
!     the relative CISE root given in the input by the 'CIRF' keyword.

! ======================================================================
! NN.14 FIXME:
!     The overlap b/w the DMRG wave of different macro iterations
!     hasn't yet been implemented. Eventually this must be fixed to
!     chose the correct root for RF calculation with DMRG-CASSCF.
!     For the time, just skip following.
! ======================================================================

!if (lRF .and. KeyCISE .and. KeyCIRF) then
!  JPCMROOT = IPCMROOT
!  IPCMROOT = IROOT(ICIRFROOT)
!  call Put_iScalar('RF CASSCF root',IPCMROOT)
!  if (JPCMROOT /= IPCMROOT) write(u6,'(1X,A,I3,A,I3)') 'RF Root has flipped from ',JPCMROOT,' to ',IPCMROOT
!else if (lRF) then
!  call Qpg_iScalar('RF CASSCF root',Exist)
!  if (.not. Exist) then
!
!    ! We are here since we are using the default values.
!
!    call Put_iScalar('RF CASSCF root',IPCMROOT)
!    call Put_iScalar('RF0CASSCF root',IPCMROOT)
!  end if
!
!  call mma_allocate(RF,nConf)
!  call Qpg_dArray('RF CASSCF Vector',Exist,mConf)
!  write(u6,*) 'Exist=',Exist
!  if (Exist .and. (mConf == nConf) .and. (iFinal /= 2) .and. ((abs(RotMax) < 1.0e-3_wp) .or. KeyCISE)) then
!    call Get_dArray('RF CASSCF Vector',RF,nConf)
!    rNorm = sqrt(DDot_(nConf,RF,1,RF,1))
!    write(u6,*) 'rNorm=',rNorm
!    JPCMROOT = IPCMROOT
!    if (rNorm > 1.0e-10_wp) then
!      call mma_allocate(Temp,nConf,Label='Temp')
!      rMax = Zero
!      jDisk = IADR15(4)
!      do i=1,lRoots
!        call DDafile(JOBIPH,2,Temp,nConf,jDisk)
!        qMax = abs(DDot_(nConf,Temp,1,RF,1))
!        write(u6,*) 'qMax=',qMax
!        if ((qMax > rMax) .and. (qMax > Half)) then
!          rMax = qMax
!          JPCMROOT = i
!        end if
!      end Do
!      call mma_deallocate(Temp)
!    end if
!  else
!    JPCMROOT = IPCMROOT
!  end if
!
!  if (JPCMROOT /= IPCMROOT) then
!    write(u6,*) ' RF Root has flipped from ',IPCMROOT,' to ',JPCMROOT
!    IPCMROOT = JPCMROOT
!    call Put_iScalar('RF CASSCF root',IPCMROOT)
!  end if
!
!  jDisk = IADR15(4)
!  do i=1,IPCMROOT-1
!    call DDafile(JOBIPH,0,rdum,nConf,jDisk)
!  end do
!  call DDafile(JOBIPH,2,RF,nConf,jDisk)
!  call Put_dArray('RF CASSCF Vector',RF,nConf)
!  call mma_deallocate(RF)
!end if

end subroutine DMRGCtl

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(DMRGCtl)

#endif
