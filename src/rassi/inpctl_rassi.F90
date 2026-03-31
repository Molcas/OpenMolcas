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

subroutine INPCTL_RASSI()

#ifdef _DMRG_
use qcmaquis_interface_cfg, only: qcmaquis_param
use qcmaquis_info, only: qcmaquis_info_init, qcm_prefixes
use qcmaquis_interface_mpssi, only: qcmaquis_mpssi_init
use rasscf_global, only: doDMRG
use rassi_data, only: NASH
use Cntrl, only: NACTE
use Definitions, only: u6
#endif
use mspt2_eigenvectors, only: init_mspt2_eigenvectors
use Symmetry_Info, only: nIrrep
use rassi_global_arrays, only: ESHFT, HAM, HDIAG, JBNUM, LROOT
use rassi_data, only: ENUC, NBASF
use Cntrl, only: BINA, HEff, IBINA, IfHDia, IFHEXT, IFShft, ISTAT, MLTPLT, MXJOB, NATO, nAtoms, NBINA, NJOB, NRNATO, NSTAT, &
                 NSTATE, RefEne
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: READ_STATES
integer(kind=iwp) :: i, JOB

! get basic info from runfile
call Get_iArray('nBas',nBasF,nIrrep)
call Get_dscalar('PotNuc',ENUC)

! Read data from the ONEINT file:
call GETCNT(NATOMS)

NSTATE = 0
! Read (and do some checking) the standard input.
call READIN_RASSI()
! if there have been no states selected at this point, we need to read
! the states later from the job files.
if (NSTATE == 0) then
  READ_STATES = .true.
  do JOB=1,NJOB
    call rdjob_nstates(JOB)
  end do
  ! store the root IDs of each state
  call mma_allocate(JBNUM,nState,Label='JBNUM')
  call mma_allocate(LROOT,nState,Label='LROOT')
  LROOT(:) = 0
  do JOB=1,NJOB
    do I=0,NSTAT(JOB)-1
      JBNUM(ISTAT(JOB)+I) = JOB
    end do
  end do
else
  READ_STATES = .false.
end if

if (NATO) then
  if (NRNATO > NSTATE) NRNATO = NSTATE
end if
if (BINA) then
  do I=1,NBINA
    if (maxval(IBINA(:,I)) > NSTATE) then
      call WarningMessage(2,'State number too high in NBINA')
      call Quit_OnUserError()
    end if
  end do
end if

#ifdef _DMRG_
!> initialize DMRG interface
if (doDMRG) then
  !> initialize only the qcm file name array (one for each job) and initialize the DMRG interface later
  call qcmaquis_info_init(njob,-1,0)
end if
#endif

!> initialize eigenvector array for mspt2 hamiltonians
call init_mspt2_eigenvectors(njob,-1,0)
! Allocate a bunch of stuff
call mma_allocate(REFENE,NSTATE,Label='RefEne')
call mma_allocate(HEFF,NSTATE,NSTATE,Label='HEff')
HEff(:,:) = Zero
if (.not. IFHEXT) then
  call mma_allocate(HAM,nState,nState,Label='HAM')
  HAM(:,:) = Zero
end if
if (.not. IFSHFT) then
  call mma_allocate(ESHFT,nSTATE,Label='ESHFT')
  ESHFT(:) = Zero
end if
if (.not. IFHDIA) call mma_Allocate(HDIAG,nState,Label='HDIAG')

! Read information on the job files and check for consistency
do JOB=1,NJOB
  call RDJOB(JOB,READ_STATES)
end do

! Number of active orbitals is taken from the first JobIph. MPS-SI cannot
! handle different active spaces per JobIph, but this is checked elsewhere
#ifdef _DMRG_
if (doDMRG) then
  !> stupid info.h defines "sum", so I cannot use the intrinsic sum function here...
  qcmaquis_param%L = 0
  do i=1,nIrrep
    qcmaquis_param%L = qcmaquis_param%L+nash(i)
  end do

  ! Initialise the new MPSSI interface
  call qcmaquis_mpssi_init(qcm_prefixes,LROOT,NSTAT(1),NJOB)

  ! Check if number of active electrons is the same for all job files
  ! Otherwise, quit on error, as Dyson orbitals are not supported yet with DMRG
  if (NJOB > 1) then
    JOB = NACTE(1)
    do i=2,NJOB
      if (NACTE(i) /= JOB) then
        call WarningMessage(2,'Number of active electrons')
        write(u6,*) ' is not the same in different JOBIPH files'
        write(u6,*) ' Dyson orbitals are not yet supported in MPSSI.'
        call Quit_OnUserError()
      end if
    end do
  end if

end if
#endif

! set orbital partitioning data
call WFNSIZES_RASSI()

! Added by Ungur Liviu on 04.11.2009
! Addition of NJOB,MSJOB and MLTPLT on RunFile.

call Put_iscalar('NJOB_SINGLE',NJOB)
call Put_iscalar('MXJOB_SINGLE',MXJOB)
call Put_iArray('MLTP_SINGLE',MLTPLT,MXJOB)

call Put_iArray('NSTAT_SINGLE',NSTAT,MXJOB)
!call Put_iArray('ISTAT_SINGLE',ISTAT,MXJOB)

! .. and print it out

! Additional input processing. Start writing report.
call INPPRC()

call mma_deallocate(REFENE)
call mma_deallocate(HEff)

end subroutine INPCTL_RASSI
