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

subroutine RASSI(IRETURN)

use Basis_Info, only: nBas
use Cntrl, only: BINA, Do_SK, DQVD, DYSEXPORT, DYSO, HOP, IFHAM, IFSO, LuExc, LuOne, LuTDM, MLTPLT, NATO, NJOB, NPROP, NSTATE, &
                 ONLY_OVERLAPS, SaveDens, SONATNSTATE, SONTOSTATES, TRACK
use Fock_util_global, only: Fake_CMO2
use frenkel_global_vars, only: doCoul, doExcitonics, eNucB, vNucB
use kVectors, only: k_Vector
use Molcas, only: MxRoot
use mspt2_eigenvectors, only: deinit_mspt2_eigenvectors
use rassi_aux, only: CMO1, CMO2, DMAB, ipglob, jDisk_TDM, Job_Index, TocM
use rassi_data, only: NBASF, NBSQ, NBST, NTDMZZ
use rassi_global_arrays, only: EIGVEC, ESHFT, HAM, HDIAG, JBNUM, LROOT, SFDYS, SODYSAMPS, SODYSAMPSI, SODYSAMPSR
use Symmetry_Info, only: nIrrep, Symmetry_Info_Free
#ifdef _HDF5_
use Dens2HDF5, only: StoreDens
use mh5, only: mh5_put_dset
use RASSIWfn, only: wfn_overlap
#endif
#ifdef _DMRG_
use qcmaquis_info, only: qcmaquis_info_deinit
use qcmaquis_interface, only: qcmaquis_interface_deinit
use rasscf_global, only: doDMRG
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: IDISK, IOPT, IRC, IRETURN, ISTATE, ISY, J, JOB, JOB1, JOB2, MPLET, nh1, NSS, NZ, NZCOUL
logical(kind=iwp) :: aux2, CLOSEONE
real(kind=wp), allocatable :: DMAT(:), DYSAMPS(:,:), ENERGY(:), OCC(:), OVLP(:,:), PROP(:,:,:), SOENE(:), TDMZZ(:), USOI(:,:), &
                              USOR(:,:), VNAT(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Prologue

IRETURN = 20

call StatusLine('RASSI: ','Starting calculation')

call GETPRINTLEVEL()

! Greetings. Default settings. Initialize data sets.
call INIT_RASSI()

CLOSEONE = .false.
IRC = -1
IOPT = 0
call OPNONE(IRC,IOPT,'ONEINT',LUONE)
if (IRC /= 0) then
  write(u6,*) 'RASSI: Error opening file'
  call ABEND()
end if
CLOSEONE = .true.

! Read and check keywords etc. from stdin. Print out.
call INPCTL_RASSI()

!SVC: prepare HDF5 wavefunction file
call CRE_RASSIWFN()

!--------  RAS wave function section --------------------------
! First, in a double loop over the states, compute any matrix
! elements required for the input RASSCF state pairs.
! Needed generalized transition density matrices are computed by
! GTDMCTL. They are written on unit LUTDM.
! Needed matrix elements are computed by PROPER.
call mma_allocate(OVLP,NSTATE,NSTATE,Label='OVLP')
call mma_allocate(DYSAMPS,NSTATE,NSTATE,Label='DYSAMPS')
call mma_allocate(EigVec,nState,nState,Label='EigVec')
call mma_allocate(ENERGY,nState,Label='Energy')
call mma_allocate(TocM,NSTATE*(NSTATE+1)/2,Label='TocM')
call mma_allocate(PROP,NSTATE,NSTATE,NPROP,LABEL='Prop')
Prop(:,:,:) = Zero
DYSAMPS(:,:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Number of basis functions
NZ = 0  ! (NBAS is already used...)
do ISY=1,nIrrep
  NZ = NZ+NBASF(ISY)
end do
if (DYSO) call mma_allocate(SFDYS,nZ,nState,nState,Label='SFDYS')

aux2 = .false.
if (doCoul) then
  call mma_allocate(eNucB,mxroot*(mxroot+1)/2,Label='eNuc')
  eNucB(:) = Zero
  NZcoul = 0  ! (NBAS is already used...)
  inquire(file='AUXRFIL2',exist=aux2)
  if (aux2) then
    call NameRun('AUXRFIL2')
  else
    call NameRun('AUXRFIL1')
  end if
  call get_iArray('nBas',nBas,nIrrep)
  NZcoul = nBas(0)
  nh1 = NZcoul*(NZcoul+1)/2
  call mma_allocate(vNucB,nh1,Label='Attr PotB')
  call Get_dArray('Nuc Potential',vNucB,nh1)
  call NameRun('#Pop')    ! switch back to old RUNFILE
end if

! Loop over jobiphs:
IDISK = 0  ! Initialize disk address for TDMs.
do JOB1=1,NJOB
  do JOB2=1,JOB1

    Fake_CMO2 = JOB1 == JOB2  ! MOs1 = MOs2  ==> Fake_CMO2=.true.

    ! Compute generalized transition density matrices, as needed:
    call GTDMCTL(PROP,JOB1,JOB2,OVLP,DYSAMPS,NZ,IDISK)
  end do
end do

#ifdef _HDF5_
call mh5_put_dset(wfn_overlap,OVLP,[NSTATE,NSTATE],[0,0])
#endif
call Put_dArray('State Overlaps',OVLP,NSTATE*NSTATE)

if (TRACK) call TRACK_STATE(OVLP)
if (TRACK .or. ONLY_OVERLAPS) then

  ! Print the overlap matrix here, since MECTL is skipped
  if (IPGLOB >= 2) then
    write(u6,*)
    write(u6,*) '     OVERLAP MATRIX FOR THE ORIGINAL STATES:'
    write(u6,*)
    do ISTATE=1,NSTATE
      write(u6,'(5(1X,F15.8))') (Ovlp(j,iState),j=1,istate)
    end do
  end if
  goto 100
end if

! Property matrix elements:
call StatusLine('RASSI: ','Computing matrix elements.')
call MECTL(PROP,OVLP,HAM,ESHFT)
!                                                                      *
!***********************************************************************
!                                                                      *
! Spin-free section

!--------  SI wave function section --------------------------
! In a second section, if Hamiltonian elements were requested,
! then also a set of secular equations are solved. This gives
! a set of non-interacting, orthonormal wave functions expressed
! as linear combinations of the input wave functions. These
! results are written out, as well as the matrix elements  of
! the eigenstates, and if requested, their density matrices
! and perhaps GTDMs.

! Hamiltonian matrix elements, eigenvectors:
if (IFHAM) then
  call StatusLine('RASSI: ','Computing Hamiltonian.')
  call EIGCTL(PROP,OVLP,DYSAMPS,HAM,EIGVEC,ENERGY)
end if

! +++ J. Creutzberg, J. Norell - 2018
! Write the spin-free Dyson orbitals to .DysOrb and .molden
! files if requested
!----------------------------------------------------------------

! Bruno Tenorio, 2020. It writes now the new Dyson norms
! See e.g. dysnorm subroutine.
if (DYSEXPORT) call WRITEDYS(DYSAMPS,SFDYS,NZ,ENERGY)
! +++

! Loop over states to compute Excitonic Couplings

if (DoExcitonics) call EXCCOUPL()

!----------------------------------------------------------------------*
! Natural orbitals, if requested:
if (NATO) then
  ! CALCULATE AND WRITE OUT NATURAL ORBITALS.
  call mma_allocate(DMAT,nBSQ,Label='DMAT')
  call mma_allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
  call mma_allocate(VNAT,nBSQ,Label='VNAT')
  call mma_allocate(OCC,nBST,Label='OCC')

  call NATORB_RASSI(DMAT,TDMZZ,VNAT,OCC,EIGVEC)
  call NATSPIN_RASSI(DMAT,TDMZZ,VNAT,OCC,EIGVEC)

  call mma_deallocate(DMAT)
  call mma_deallocate(TDMZZ)
  call mma_deallocate(VNAT)
  call mma_deallocate(OCC)
end if
! Bi-natural orbitals, if requested:
if (BINA) call BINAT()

if (.not. IFHAM) goto 100

!                                                                      *
!***********************************************************************
!                                                                      *
!  Spin_Orbit section

!-------- Spin-Orbit calculations   --------------------------
! In this third section, if spin-orbit coupling parameters were
! computed by the previous sections, then additionally the
! spin-orbit eigenfunctions and levels are computed.

! Nr of spin states and division of loops:
NSS = 0
do ISTATE=1,NSTATE
  JOB = JBNUM(ISTATE)
  MPLET = MLTPLT(JOB)
  NSS = NSS+MPLET
end do

call mma_allocate(USOR,NSS,NSS,Label='USOR')
call unitmat(USOR,NSS)
call mma_allocate(USOI,NSS,NSS,Label='USOI')
USOI(:,:) = Zero
call mma_allocate(SOENE,nSS,Label='SOENE')
SOENE(:) = Zero

if (IFSO) then
  call StatusLine('RASSI: ','Computing SO Hamiltonian.')
  call SOEIG(PROP,USOR,USOI,SOENE,NSS,ENERGY)
end if

#ifdef _HDF5_
! Store TDMs in HDF5
call StoreDens(EigVec)
#endif

! +++ J. Norell - 2018
! Make the SO Dyson orbitals and amplitudes from the SF ones

if (DYSO) then

  if (IFSO) then
    call mma_allocate(SODYSAMPS,NSS,NSS,Label='SODYSAMPS')
    call mma_allocate(SODYSAMPSR,NSS,NSS,Label='SODYSAMPSR')
    call mma_allocate(SODYSAMPSI,NSS,NSS,Label='SODYSAMPSI')

    call SODYSORB(NSS,USOR,USOI,DYSAMPS,NZ,SOENE)
  end if

  call mma_deallocate(SFDYS)
end if

! +++

call PRPROP(PROP,USOR,USOI,SOENE,NSS,OVLP,ENERGY,JBNUM,EigVec)

! Plot SO-Natural Orbitals if requested
! Will also handle mixing of states (sodiag)
if (SONATNSTATE > 0) call DO_SONATORB(NSS,USOR,USOI)
! Plot SO-Natural Transition Orbitals if requested
if (SONTOSTATES > 0) call DO_SONTO(NSS,USOR,USOI)

call mma_deallocate(USOR)
call mma_deallocate(USOI)
call mma_deallocate(SOENE)
if (DYSO .and. IFSO) then
  call mma_deallocate(SODYSAMPS)
  call mma_deallocate(SODYSAMPSR)
  call mma_deallocate(SODYSAMPSI)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Trajectory Surface Hopping

! Turns on the procedure if the Keyword HOP was specified.

if (HOP) then
  call StatusLine('RASSI: ','Trajectory Surface Hopping')
  call TSHinit(ENERGY(:))
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! CEH April 2015 DQV diabatization scheme
! This passes the PROP matrix into the DQV diabatization subroutine
! The subroutine will compute a transformation matrix, which is used
! to compute diabats.

! The user has to compute x, y, z, xx, yy, zz, 1/r in this order.

if (DQVD) then
  call StatusLine('RASSI: ','DQV Diabatization')
  call DQVDiabat(PROP,HAM)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!     EPILOGUE                                                         *
!                                                                      *
!***********************************************************************
!                                                                      *
100 continue

if (DoCoul) then
  if (.not. aux2) then
    call NameRun('AUXRFIL1')
    call Put_dArray('<rhoB|VnucA>',eNucB,mxroot*(mxroot+1)/2)
    call Cho_X_Final(irc)
    call NameRun('#Pop') ! switch back to old RUNFILE
    call mma_deallocate(VNucB)
    call mma_deallocate(eNucB)
  else
    call NameRun('AUXRFIL1')
    call Cho_X_Final(irc)
    call NameRun('#Pop')
    call mma_deallocate(VNucB)
    call mma_deallocate(eNucB)
  end if
end if

call mma_deallocate(Ovlp)
call mma_deallocate(DYSAMPS)
call mma_deallocate(HAM)
call mma_deallocate(EigVec)
call mma_deallocate(Energy)
call mma_deallocate(ESHFT)
call mma_deallocate(HDIAG)
call mma_deallocate(jDisk_TDM)
call mma_deallocate(JBNUM)
call mma_deallocate(LROOT)
call mma_deallocate(TocM)
call mma_deallocate(Prop)

if (Do_SK) call mma_deallocate(k_Vector)

if (CLOSEONE) then
  IRC = -1
  call CLSONE(IRC,0)
  if (IRC /= 0) then
    write(u6,*) 'RASSI: Error opening file'
    call ABEND()
  end if
end if

#ifdef _DMRG_
!> finalize MPS-SI interface
if (doDMRG) then
  call qcmaquis_interface_deinit()
  call qcmaquis_info_deinit()
end if
#endif
!> free memory (if allocated at all - currently only for QD-NEVPT2 as ref wfn)
call deinit_mspt2_eigenvectors()
!                                                                      *
!***********************************************************************
!                                                                      *
! Close dafiles.

if (SaveDens) then
  call DaClos(LuTDM)
  call mma_deallocate(JOB_INDEX,safe='*')
  call mma_deallocate(CMO1,safe='*')
  call mma_deallocate(CMO2,safe='*')
  call mma_deallocate(DMAB,safe='*')
end if
call DaClos(LuExc)
call Symmetry_Info_Free()
!                                                                      *
!***********************************************************************
!                                                                      *
! PRINT I/O STATISTICS:
call FASTIO('STATUS')

call StatusLine('RASSI: ','Finished.')
IRETURN = 0

end subroutine RASSI
