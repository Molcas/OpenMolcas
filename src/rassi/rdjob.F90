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

subroutine RDJOB(JOB,READ_STATES)

use Cntrl, only: bNAME, HAVE_DIAG, HAVE_HEFF, HEAD1, HEff, IDCMO, IFEJOB, IFHEFF, IROOT1, IRREP, ISTAT, iTOC15, JBNAME, LROT1, &
                 LSYM1, LuIPH, MLTPLT, MPLET1, NACTE, NACTE1, NASH1, NBAS1, NCONF, NCONF1, NDEL1, NELE3, NELE31, NFRO1, NHOL11, &
                 NHOLE1, NISH1, NROOT1, NRS11, NRS21, NRS31, NSTAT, NSTATE, NSYM1, RASTYP, RefEne, TITLE1
use gugx, only: LEVEL
use Molcas, only: LenIn, MxOrb, MxRoot, MxSym, MxLev
use rasdef, only: NRS1, NRS2, NRS3
use RASDim, only: MxIter, MxTit
use rassi_aux, only: ipglob
use rassi_data, only: WFTYPE, NASH, NSSH, NDEL, NOSH, NASH, NISH, NFRO, NBASF, NDEL, NFRO, NISH
use rassi_global_arrays, only: JBNUM, LROOT
use Symmetry_Info, only: nIrrep
#ifdef _HDF5_
use Cntrl, only: NDET, NROOTS, QDPT2SC
use mh5, only: mh5_close_file, mh5_exists_attr, mh5_exists_dset, mh5_fetch_attr, mh5_fetch_dset, mh5_is_hdf5, mh5_open_file_r
#endif
#ifdef _DMRG_
use qcmaquis_info, only: qcm_group_names, qcm_prefixes, qcmaquis_info_init
use rasscf_global, only: doDMRG
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: JOB
logical(kind=iwp) :: READ_STATES
integer(kind=iwp) :: I, IAD, IAD15, IDISK, IDUM(1), IERR, IPT2, ISNUM, ISTATE, ISY, IT, J, JSNUM, JSTATE, NEJOB, NHEFF, NIS, NIS1, &
                     NMAYBE, NROOT0, NTIT1
real(kind=wp) :: AEMAX, E, ENUCDUMMY, HIJ
logical(kind=iwp) :: ISZERO
real(kind=wp), allocatable :: EJOB(:), EREAD(:), H_Eff(:), Weight(:)
#ifdef _HDF5_
integer(kind=iwp) :: ref_iSpin, ref_nactel, ref_nash(mxSym), ref_nBas(mxSym), ref_nconf, ref_ndel(mxSym), ref_ndet, ref_nelec3, &
                     ref_nfro(mxSym), ref_nhole1, ref_nish(mxSym), ref_nroots, ref_nrs1(mxSym), ref_nrs2(mxSym), ref_nrs3(mxSym), &
                     ref_nssh(mxSym), ref_nstates, ref_nSym, ref_stSym, refwfn_id
character(len=21) :: pt2_e_string
character(len=16) :: molcas_module
character(len=8) :: heff_string
integer(kind=iwp), allocatable :: ref_rootid(:), root2state(:)
real(kind=wp), allocatable :: ref_energies(:), ref_Heff(:,:)
character, allocatable :: typestring(:)
#endif
#ifdef _DMRG_
integer(kind=iwp) :: idx
character(len=300) :: WorkDir
#endif

#ifdef _HDF5_
if (mh5_is_hdf5(jbname(job))) then
  !*********************************************************************
  !
  ! For HDF5 formatted job files
  !
  !*********************************************************************

  if (IPGLOB >= 2) then
    if (JOB == 1) then
      write(u6,*)
      write(u6,'(6X,A)') repeat('*',80)
      write(u6,'(6X,A1,78X,A1)') '*','*'
      write(u6,'(6X,A1,24X,A,24X,A1)') '*','     General data section     ','*'
      write(u6,'(6X,A1,78X,A1)') '*','*'
      write(u6,'(6X,A)') repeat('*',80)
    end if
    write(u6,*)
    write(u6,*) '  Specific data for HDF5 file ',trim(JBNAME(JOB))
    write(u6,*) '  -------------------------------------'
  end if

  refwfn_id = mh5_open_file_r(jbname(job))

  call mh5_fetch_attr(refwfn_id,'MOLCAS_MODULE',molcas_module)
  call mh5_fetch_attr(refwfn_id,'SPINMULT',ref_iSpin)
  call mh5_fetch_attr(refwfn_id,'NSYM',ref_nSym)
  call mh5_fetch_attr(refwfn_id,'LSYM',ref_stSym)
  call mh5_fetch_attr(refwfn_id,'NBAS',ref_nBas)

  call mh5_fetch_attr(refwfn_id,'NACTEL',ref_nactel)
  call mh5_fetch_attr(refwfn_id,'NHOLE1',ref_nhole1)
  call mh5_fetch_attr(refwfn_id,'NELEC3',ref_nelec3)
  call mh5_fetch_attr(refwfn_id,'NCONF',ref_nconf)
  call mh5_fetch_attr(refwfn_id,'NSTATES',ref_nstates)
  if (mh5_exists_attr(refwfn_id,'NROOTS')) then
    call mh5_fetch_attr(refwfn_id,'NROOTS',ref_nroots)
  else
    ref_nroots = ref_nstates
  end if
  ! NDET array is read only from HDF5, the number is not in JOBIPH
  if (mh5_exists_attr(refwfn_id,'NDET')) then
    call mh5_fetch_attr(refwfn_id,'NDET',ref_ndet)
  else
    ! to avoid runtime error
    ref_ndet = 1
  end if

  call mma_allocate(typestring,sum(ref_nbas(1:ref_nsym)))
  call mh5_fetch_dset(refwfn_id,'MO_TYPEINDICES',typestring)
  call tpstr2orb(ref_nsym,ref_nbas,typestring,ref_nfro,ref_nish,ref_nrs1,ref_nrs2,ref_nrs3,ref_nssh,ref_ndel)
  ref_nash(1:nIrrep) = ref_nrs1(1:nIrrep)+ref_nrs2(1:nIrrep)+ref_nrs3(1:nIrrep)
  call mma_deallocate(typestring)

# ifdef _DMRG_
  if ((.not. mh5_exists_dset(refwfn_id,'CI_VECTORS')) .and. (.not. doDMRG)) then
    ! Leon: TODO: This must be also extended for other DMRG interfaces
    ! than QCMaquis
# else
  if (.not. mh5_exists_dset(refwfn_id,'CI_VECTORS')) then
# endif
    write(u6,'(1X,A)') 'The HDF5 file does not contain CI vectors,'
    write(u6,'(1X,A)') 'make sure it was created by rasscf/caspt2.'
    call AbEnd()
  end if
  if (.not. mh5_exists_dset(refwfn_id,'MO_VECTORS')) then
    write(u6,'(1X,A)') 'The HDF5 file does not contain MO vectors,'
    write(u6,'(1X,A)') 'make sure it was created by rasscf/caspt2/nevpt2.'
    call AbEnd()
  end if

  !call mh5_fetch_attr (refwfn_id,'L2ACT', L2ACT)
  call mh5_fetch_attr(refwfn_id,'A2LEV',LEVEL)

  call mma_allocate(ref_rootid,ref_nstates)
  call mh5_fetch_attr(refwfn_id,'STATE_ROOTID',ref_rootid)
  call mma_allocate(root2state,MxRoot,Label='root2state')
  call iCopy(MxRoot,[0],0,root2state,1)
  if (mh5_exists_attr(refwfn_id,'ROOT2STATE')) then
    call mh5_fetch_attr(refwfn_id,'ROOT2STATE',root2state)
  else
    do i=1,ref_nroots
      root2state(i) = i
    end do
  end if
  if (read_states) then
    !  Do not update the state number here, because it's already read in
    !  rdjob_nstates()
    !        NSTAT(JOB)=ref_nstates
    !        NSTATE=NSTATE+ref_nstates
    ! store the root IDs of each state
    do I=0,NSTAT(JOB)-1
      LROOT(ISTAT(JOB)+I) = ref_rootid(I+1)
      JBNUM(ISTAT(JOB)+I) = JOB
    end do
  end if
  LROT1 = ref_nroots
  do I=0,NSTAT(JOB)-1
    NROOT0 = root2state(LROOT(ISTAT(JOB)+I))
    if ((NROOT0 <= 0) .or. (NROOT0 > LROT1)) goto 9002
  end do

  if (qdpt2sc .and. (molcas_module(1:6) == 'NEVPT2')) then
    heff_string = 'H_EFF_SC'
    pt2_e_string = 'STATE_PT2_ENERGIES_SC'
  else
    heff_string = 'H_EFF'
    pt2_e_string = 'STATE_PT2_ENERGIES'
  end if

  if (mh5_exists_dset(refwfn_id,heff_string)) then
    ! read the ms-caspt2/qd-nevpt2 effective hamiltonian if it is available
    call mma_allocate(ref_Heff,ref_nstates,ref_nstates)
    call mh5_fetch_dset(refwfn_id,heff_string,ref_Heff)
    HAVE_HEFF = .true.
    ! with ejob, only read diagonal
    if (ifejob) then
      HAVE_DIAG = .true.
      !call WarningMessage(0,'Effective Hamiltonian found in  reference file, but "EJOB" was requested: off-diagonal elements '// &
      !                    'will be ignored!')
      do I=1,NSTAT(JOB)
        ISTATE = ISTAT(JOB)-1+I
        ISNUM = root2state(LROOT(ISTATE))
        REFENE(istate) = ref_Heff(ISNUM,ISNUM)
      end do
    else
      write(u6,'(2x,a)') ' Effective Hamiltonian from MRPT2 in action'
      write(u6,'(2x,a)') ' ------------------------------------------'
      do I=1,NSTAT(JOB)
        ISTATE = ISTAT(JOB)-1+I
        ISNUM = root2state(LROOT(ISTATE))
        do J=1,NSTAT(JOB)
          JSTATE = ISTAT(JOB)-1+J
          JSNUM = root2state(LROOT(JSTATE))
          HEff(jState,iState) = ref_Heff(ISNUM,JSNUM)
        end do
      end do
    end if
    call mma_deallocate(ref_Heff)
  else if (mh5_exists_dset(refwfn_id,pt2_e_string)) then
    ! read the caspt2/qdnevpt2 reference energies if available
    HAVE_DIAG = .true.
    call mma_allocate(ref_energies,ref_nstates)
    call mh5_fetch_dset(refwfn_id,pt2_e_string,ref_energies)
    do I=1,NSTAT(JOB)
      ISTATE = ISTAT(JOB)-1+I
      ISNUM = root2state(LROOT(ISTATE))
      REFENE(istate) = ref_energies(ISNUM)
      ! put the energies on the Heff diagonal too, just in case
      HEff(iState,iState) = ref_energies(ISNUM)
    end do
    call mma_deallocate(ref_energies)
  else if (mh5_exists_dset(refwfn_id,'ROOT_ENERGIES')) then
    ! read rasscf energies
    HAVE_DIAG = .true.
    call mma_allocate(ref_energies,ref_nroots)
    call mh5_fetch_dset(refwfn_id,'ROOT_ENERGIES',ref_energies)
    do I=1,NSTAT(JOB)
      ISTATE = ISTAT(JOB)-1+I
      ISNUM = root2state(LROOT(ISTATE))
      REFENE(istate) = ref_energies(ISNUM)
      ! put the energies on the Heff diagonal too, just in case
      HEff(iState,iState) = ref_energies(ISNUM)
    end do
    call mma_deallocate(ref_energies)
  end if
  call mma_deallocate(root2state)

# ifdef _DMRG_
  call getenvf('WorkDir',WorkDir)
  ! Leon 5/12/2016: Fetch QCMaquis checkpoint names if requested
  if (doDMRG) then
    if (mh5_exists_dset(refwfn_id,'QCMAQUIS_CHECKPOINT')) then
      write(u6,*) '  QCMaquis checkpoint files:'
      write(u6,*) '  --------------------------'
      write(u6,*) '  State   Checkpoint file'

      !> allocate space for the file name strings of job JOB
      call qcmaquis_info_init(job,nstat(job),1)

      do I=1,NSTAT(JOB)
        ISTATE = ISTAT(JOB)-1+I
        call mh5_fetch_dset(refwfn_id,'QCMAQUIS_CHECKPOINT',qcm_group_names(job)%states(i:i),[1],[LROOT(ISTATE)-1])
        write(u6,'(5X,I3,3X,A)') ISTATE,trim(qcm_group_names(job)%states(i))
      end do
      write(u6,*) '  --------------------------'
      ! save QCMaquis prefix
      ! by cutting off the last '.checkpoint_state.X.h5'
      ! and adding the full path
      if (size(qcm_group_names(job)%states) > 0) then
        idx = index(qcm_group_names(job)%states(1),'.checkpoint_state.')
        if (idx > 0) then
          qcm_prefixes(job) = trim(WorkDir)//'/'//trim(qcm_group_names(job)%states(1)(1:idx-1))
        else
          call WarningMessage(2,'Faulty QCMaquis checkpoint name')
          write(u6,*) 'Must contain "checkpoint_state"'
          call Abend()
        end if
      end if
    else
      call WarningMessage(2,'QCMaquis checkpoint names not found on HDF5 files. Make sure you created them with the MOLCAS '// &
                          'version which supports them')
      call Quit_OnUserError()
    end if
  end if
# endif
  if (ref_nsym /= nIrrep) then
    call WarningMessage(2,'NSYM not consistent with RunFile')
    call Quit_OnUserError()
  end if
  do i=1,nIrrep
    if (ref_nbas(i) /= nbasf(i)) then
      call WarningMessage(2,'NBAS not consistent with RunFile')
      call Quit_OnUserError()
    end if
  end do

  NACTE(JOB) = ref_nactel
  NHOLE1(JOB) = ref_nhole1
  NELE3(JOB) = ref_nelec3
  MLTPLT(JOB) = ref_iSpin
  IRREP(JOB) = ref_stSym
  NCONF(JOB) = ref_nConf
  NROOTS(JOB) = ref_nroots
  ! in singlet case the number of determinants is doubled in rassi
  ! compare to the rasscf routine, storing here due to rassi procedure
  if (mltplt(JOB) == 1) then
    nDet(JOB) = 2*ref_ndet-1
  else
    nDet(JOB) = ref_ndet
  end if

  if (job == 1) then
    ! first wavefunction file, set global variables
    do I=1,nIrrep
      NFRO(I) = 0
      NISH(I) = ref_nfro(I)+ref_nish(I)
      NASH(I) = ref_nash(I)
      NRS1(I) = ref_NRS1(I)
      NRS2(I) = ref_NRS2(I)
      NRS3(I) = ref_NRS3(I)
      NOSH(I) = NISH(I)+NASH(I)
      NDEL(I) = 0
      NSSH(I) = NBASF(I)-NFRO(I)-NISH(I)-NASH(I)-NDEL(I)
    end do
  else
    ! subsequent wavefunction file, check against global variables
    if ((ref_nhole1 /= nhole1(1)) .or. (ref_nelec3 /= nele3(1))) then
      call WarningMessage(2,'inconsistent RAS holes/electrons')
      call Quit_OnUserError()
    end if
    do i=1,nIrrep
      if ((ref_nfro(i)+ref_nish(i) /= nish(i)) .or. (ref_nash(i) /= nash(i)) .or. (ref_nrs1(i) /= nrs1(i)) .or. &
          (ref_nrs2(i) /= nrs2(i)) .or. (ref_nrs3(i) /= nrs3(i))) then
        call WarningMessage(2,'inconsistent orbital partitioning')
        call Quit_OnUserError()
      end if
    end do
  end if

  WFTYPE = 'GENERAL'
  if (ref_nactel == 2*sum(NASH(1:nIrrep))) WFTYPE = 'CLOSED'
  if (ref_nactel == 0) WFTYPE = 'EMPTY'
  RASTYP(JOB) = WFTYPE

  if (IPGLOB >= 2) then
    write(u6,'(A,I9)') '  STATE IRREP:        ',IRREP(JOB)
    write(u6,'(A,I9)') '  SPIN MULTIPLICITY:  ',MLTPLT(JOB)
    write(u6,'(A,I9)') '  ACTIVE ELECTRONS:   ',NACTE(JOB)
#   ifdef _DMRG_
    if (.not. doDMRG) then
#   endif
      write(u6,'(A,I9)') '  MAX RAS1 HOLES:     ',NHOLE1(JOB)
      write(u6,'(A,I9)') '  MAX RAS3 ELECTRONS: ',NELE3(JOB)
      write(u6,'(A,I9)') '  NR OF CONFIG:       ',NCONF(JOB)
#   ifdef _DMRG_
    end if
#   endif
  end if
  if (IPGLOB >= 3) write(u6,*) '  Wave function type WFTYPE=',WFTYPE

  call mma_deallocate(ref_rootid)
  call mh5_close_file(refwfn_id)

else
#endif
  !*********************************************************************
  !
  ! For JOBIPH/JOBMIX formatted job files
  !
  !*********************************************************************
# ifdef _DMRG_
  if (doDMRG) then
    call WarningMessage(3,'QCMaquis requires checkpoint names from JOBxxx. This works only with HDF5 JobIph files. Please make '// &
                        'sure you use a .h5 file as JOBxxx.')
    call abend()
  end if
# endif
  if (IPGLOB >= 2) then
    if (JOB == 1) then
      write(u6,*)
      write(u6,'(6X,A)') repeat('*',80)
      write(u6,'(6X,A1,78X,A1)') '*','*'
      write(u6,'(6X,A1,24X,A,24X,A1)') '*','     General data section     ','*'
      write(u6,'(6X,A1,78X,A1)') '*','*'
      write(u6,'(6X,A)') repeat('*',80)
    end if
    write(u6,*)
    write(u6,*) '  Specific data for JOBIPH file ',trim(JBNAME(JOB))
    write(u6,*) '  -------------------------------------'
  end if
  ! Open JOBIPH file:
  call DANAME(LUIPH,JBNAME(JOB))
  ! READ TABLE OF CONTENTS ON THIS JOBIPH FILE:
  IAD = 0
  call IDAFILE(LUIPH,2,ITOC15,30,IAD)
  ! SCATTER-READ VARIOUS DATA:
  ! PAM Mar2014: Note that ENUC1 (=POTNUC) replaced by dummy placeholder
  IAD = ITOC15(1)
  call mma_allocate(Weight,MxRoot,Label='Weight')
  call WR_RASSCF_Info(LUIPH,2,IAD,NACTE1,MPLET1,NSYM1,LSYM1,NFRO1,NISH1,NASH1,NDEL1,NBAS1,mxSym,bNAME,(LenIn+8)*mxOrb,NCONF1, &
                      HEAD1,2*72,TITLE1,4*mxTit*18,ENUCDUMMY,LROT1,NROOT1,IROOT1,mxRoot,NRS11,NRS21,NRS31,NHOL11,NELE31,IPT2,Weight)
  call mma_deallocate(Weight)
  ! Response field contribution to zero-electron energies
  ! is added in GETH1.
  if (READ_STATES) then
    ! Do not update the state number here, because it's already read in
    ! rdjob_nstates()
    !ISTAT(JOB) = NSTATE+1
    !NSTAT(JOB) = NROOT1
    !NSTATE = NSTATE+NROOT1
    ! store the root IDs of each state

    ! If unset yet, set now
    if (LROOT(ISTAT(JOB)) == 0) then
      do I=0,NSTAT(JOB)-1
        LROOT(ISTAT(JOB)+I) = IROOT1(I+1)
        JBNUM(ISTAT(JOB)+I) = JOB
      end do
    end if
  end if
  do I=0,NSTAT(JOB)-1
    NROOT0 = LROOT(ISTAT(JOB)+I)
    if (NROOT0 > LROT1) goto 9002
  end do

  ! First read energies, which may be used in any case
  NEJOB = MXROOT*MXITER
  call mma_allocate(EJOB,NEJOB,Label='EJOB')
  IAD = ITOC15(6)
  call DDAFILE(LUIPH,2,EJOB,NEJOB,IAD)
  ! Note that there is no info on nr of iterations
  ! so we cannot know what energies to pick...
  ! Let us make a guess: The correct set of energy values in the
  ! table of energies/iteration is the last one with not all zeroes.
  NMAYBE = 0
  do IT=1,MXITER
    AEMAX = Zero
    do I=1,MXROOT
      E = EJOB(MXROOT*(IT-1)+I)
      AEMAX = max(AEMAX,abs(E))
    end do
    if (abs(AEMAX) <= 1.0e-12_wp) exit
    NMAYBE = IT
  end do
  call mma_allocate(EREAD,NSTATE,Label='EREAD')
  do I=1,NSTAT(JOB)
    ISTATE = ISTAT(JOB)-1+I
#   ifdef _DMRG_
    if (doDMRG) then
      E = EJOB(LROOT(ISTATE)-ISTAT(JOB)+1+MXROOT*(NMAYBE-1))
    else
#   endif
      E = EJOB(LROOT(ISTATE)+MXROOT*(NMAYBE-1))
#   ifdef _DMRG_
    end if
#   endif
    EREAD(istate) = E
  end do
  call mma_deallocate(EJOB)

  ! Using energy data from JobIph?
  if (IFEJOB) then
    if (ITOC15(15) == -1) HAVE_HEFF = .true.
    if (NMAYBE == 0) then
      write(u6,*) ' Sorry. Keyword ''EJOB'' has been used'
      write(u6,*) ' but there are no energies available on'
      write(u6,*) ' the JOBIPH file nr',JOB
      call ABEND()
    end if
    HAVE_DIAG = .true.
    ! Put the energies into diagonal of Hamiltonian:
    do I=1,NSTAT(JOB)
      ISTATE = ISTAT(JOB)-1+I
      REFENE(istate) = EREAD(istate)
    end do
  end if

  ! Using effective Hamiltonian from JobIph file?
  if (IFHEFF) then
    if (ITOC15(15) /= -1) then
      write(u6,*) 'RDJOB Error: HEFF not found on JOBIPH.'
      write(u6,*) 'The HEFF keyword was used, but the JOBIPH file'
      write(u6,*) 'uses an old layout where this data field is'
      write(u6,*) 'not present. Recompute JOBIPH file, or put'
      write(u6,*) 'effective Hamiltonian in input after keyword'
      write(u6,*) 'HEXT. Program stops here.'
      call ABEND()
    end if
    HAVE_HEFF = .true.
    NHEFF = LROT1**2
    call mma_allocate(H_EFF,NHEFF,Label='H_Eff')
    IAD15 = ITOC15(17)
    call DDAFILE(LUIPH,2,H_EFF,NHEFF,IAD15)
    ! If both EJOB and HEFF are given, read only the diagonal
    if (IFEJOB) then
      HAVE_DIAG = .true.
      do I=1,NSTAT(JOB)
        ISTATE = ISTAT(JOB)-1+I
        ISNUM = LROOT(ISTATE)
        HIJ = H_EFF(ISNUM+LROT1*(ISNUM-1))
        REFENE(istate) = HIJ
        HEff(iState,iState) = HIJ
      end do
    else
      do I=1,NSTAT(JOB)
        ISTATE = ISTAT(JOB)-1+I
        ISNUM = LROOT(ISTATE)
        ISZERO = .true.
        do J=1,NSTAT(JOB)
          JSTATE = ISTAT(JOB)-1+J
          JSNUM = LROOT(JSTATE)
          HIJ = H_EFF(ISNUM+LROT1*(JSNUM-1))
          HEff(jState,iState) = HIJ
          if (I == J) REFENE(istate) = HIJ
          if (abs(HIJ) > Zero) ISZERO = .false.
        end do
        if (ISZERO) then
          Heff(iState,iState) = EREAD(istate)
        end if
      end do
    end if
    call mma_deallocate(H_Eff)
  end if
  call mma_deallocate(EREAD)
  ! Read the level to orbital translations
  IDISK = ITOC15(18)
  call IDAFILE(LUIPH,0,IDUM,MXLEV,IDISK) ! L2ACT
  call IDAFILE(LUIPH,2,LEVEL,MXLEV,IDISK)
  ! Close JobIph file
  call DACLOS(LUIPH)

  ! The RASSCF program is not certain to give consistent data. For
  ! pure CASSCF cases, it may not bother to set the NRS1..NRS3 arrays.
  ! Check and repair:
  if (NHOL11+NELE31 == 0) then
    IERR = 0
    do I=1,NSYM1
      if (NRS11(I) /= 0) IERR = 1
      if (NRS21(I) /= NASH1(I)) IERR = 1
      if (NRS31(I) /= 0) IERR = 1
    end do
    if (IERR == 1) then
      write(u6,*)
      write(u6,*) ' (NOTE: The nr of RAS1, RAS2 and RAS3 orbitals as recorded on the JOBIPH file do not match'
      write(u6,*) ' the number of active orbitals. But this is a pure CASSCF case. Maybe the RASSCF programmer did'
      write(u6,*) ' not bother with the RAS1..RAS3 arrays in that case. RASSI will reset these arrays as needed.)'
      write(u6,*)
      do I=1,NSYM1
        NRS11(I) = 0
        NRS21(I) = NASH1(I)
        NRS31(I) = 0
      end do
    end if
  end if

  if (JOB == 1) then
    ! FIRST JOB FILE. TRANSFER DATA TO COMMON:
    nIrrep = NSYM1
    do I=1,nIrrep
      NFRO(I) = 0
      NISH(I) = NFRO1(I)+NISH1(I)
      NASH(I) = NASH1(I)
      NRS1(I) = NRS11(I)
      NRS2(I) = NRS21(I)
      NRS3(I) = NRS31(I)
      NOSH(I) = NISH(I)+NASH(I)
      NDEL(I) = 0
      NBASF(I) = NBAS1(I)
      NSSH(I) = NBASF(I)-NFRO(I)-NISH(I)-NASH(I)-NDEL(I)
    end do
  else
    ! THIS IS NOT THE FIRST JOBIPH.
    ! CHECK THAT DATA IS CONSISTENT WITH EARLIER:
    if (NSYM1 /= nIrrep) goto 9001
    if (NHOL11 /= NHOLE1(JOB-1)) goto 9003
    if (NELE31 /= NELE3(JOB-1)) goto 9003
    do ISY=1,NSYM1
      NIS1 = NISH1(ISY)+NFRO1(ISY)
      NIS = NISH(ISY)
      if (NIS1 /= NIS) goto 9004
      if (NRS11(ISY) /= NRS1(ISY)) goto 9005
      if (NRS21(ISY) /= NRS2(ISY)) goto 9005
      if (NRS31(ISY) /= NRS3(ISY)) goto 9005
      if (NBAS1(ISY) /= NBASF(ISY)) goto 9006
    end do
  end if

  ! DATA PARTICULAR TO THIS JOBIPH:
  if (IPGLOB >= 2) then
    write(u6,*)
    write(u6,*) '  Header from SEWARD:'
    write(u6,'(7X,36A2)') (HEAD1(I),I=1,36)
    write(u6,'(7X,36A2)') (HEAD1(I),I=37,72)
    ! NOTE: AT PRESENT, JOBIPH FILE GIVES NO INFORMATION ON THE
    ! AMOUNT OF TITLE LINES.
    NTIT1 = 1
    write(u6,*)
    write(u6,*) '  CASSCF title (first line only):'
    write(u6,'(7X,18A4)') ((TITLE1(I,J),I=1,18),J=1,NTIT1)
    write(u6,*)
    write(u6,'(A,I9)') '  STATE IRREP:        ',LSYM1
    write(u6,'(A,I9)') '  SPIN MULTIPLICITY:  ',MPLET1
    write(u6,'(A,I9)') '  ACTIVE ELECTRONS:   ',NACTE1
    write(u6,'(A,I9)') '  MAX RAS1 HOLES:     ',NHOL11
    write(u6,'(A,I9)') '  MAX RAS3 ELECTRONS: ',NELE31
    write(u6,'(A,I9)') '  NR OF CONFIG:       ',NCONF1
  end if
  WFTYPE = 'GENERAL'
  !if (MPLET1 == (SUM(NASH(1:nIrrep))+1)) WFTYPE = 'HISPIN'
  ! Note: the HISPIN case may be buggy and is not used presently.
  if (MPLET1 == (sum(NASH(1:nIrrep))+1)) then
    write(u6,*) ' This wave function is of HISPIN type.'
    write(u6,*) ' However, the special handling for that case'
    write(u6,*) ' is suspected to be buggy. So the variable'
    write(u6,*) ' WFTYPE is set to GENERAL.'
  end if
  if (NACTE1 == 2*sum(NASH(1:nIrrep))) WFTYPE = 'CLOSED'
  if (NACTE1 == 0) WFTYPE = 'EMPTY'
  RASTYP(JOB) = WFTYPE
  if (IPGLOB >= 3) write(u6,*) '  Wave function type WFTYPE=',WFTYPE
  NACTE(JOB) = NACTE1
  MLTPLT(JOB) = MPLET1
  IRREP(JOB) = LSYM1
  NHOLE1(JOB) = NHOL11
  NELE3(JOB) = NELE31
  NCONF(JOB) = NCONF1
  ! Where is the CMO data set stored?
  IDCMO(JOB) = ITOC15(2)
  if (IPT2 /= 0) IDCMO(JOB) = ITOC15(9)

#ifdef _HDF5_
end if
#endif

call XFLUSH(u6)

return

!***********************************************************************
!
! Error exits
!
!***********************************************************************
9001 write(u6,*) ' SYMMETRY GROUPS MUST BE EQUAL.'
write(u6,*) ' NSYM1:',NSYM1,'NSYM :',nIrrep
goto 9010
9002 write(u6,*) ' ROOT NOT AVAILABLE.'
write(u6,*) '             REQUESTED ROOT:',NROOT0
write(u6,*) '  MAXIMUM ROOT IN THIS FILE:',LROT1
goto 9010
9003 write(u6,*) ' RAS SPECIFICATIONS DIFFER.'
write(u6,*) '     THIS STATE: MAX NR OF RAS-1 HOLES:',NHOL11
write(u6,*) '             MAX NR OF RAS-3 ELECTRONS:',NELE31
write(u6,*) ' PREVIOUS STATE: MAX NR OF RAS-1 HOLES:',NHOLE1(JOB-1)
write(u6,*) '             MAX NR OF RAS-3 ELECTRONS:',NELE3(JOB-1)
goto 9010
9004 write(u6,*) ' NR. OF (FROZEN+INACTIVE) ORBITALS DIFFER.'
write(u6,'(A,8I4)') ' THIS STATE:',(NFRO1(I)+NISH1(I),I=1,NSYM1)
write(u6,'(A,8I4)') '   PREVIOUS:',(NISH(I),I=1,nIrrep)
goto 9010
9005 write(u6,*) ' NR. OF ACTIVE ORBITALS DIFFER.'
write(u6,'(A,8I4)') ' THIS STATE, ACTIVE:',(NASH1(I),I=1,NSYM1)
write(u6,'(A,8I4)') '               RAS1:',(NRS11(I),I=1,NSYM1)
write(u6,'(A,8I4)') '               RAS2:',(NRS21(I),I=1,NSYM1)
write(u6,'(A,8I4)') '               RAS3:',(NRS31(I),I=1,NSYM1)
write(u6,'(A,8I4)') '   PREVIOUS, ACTIVE:',(NASH(I),I=1,nIrrep)
write(u6,'(A,8I4)') '               RAS1:',(NRS1(I),I=1,NSYM1)
write(u6,'(A,8I4)') '               RAS2:',(NRS2(I),I=1,NSYM1)
write(u6,'(A,8I4)') '               RAS3:',(NRS3(I),I=1,NSYM1)
goto 9010
9006 write(u6,*) ' NR. OF BASIS FUNCTION DIFFER.'
write(u6,'(A,8I4)') ' THIS JOBIPH:',(NBAS1(I),I=1,NSYM1)
write(u6,'(A,8I4)') '    PREVIOUS:',(NBASF(I),I=1,nIrrep)
9010 continue
write(u6,*) ' DATA IN JOBIPH FILE NAMED ',trim(JBNAME(JOB)),' WERE'
write(u6,*) ' INCONSISTENT WITH EARLIER DATA. PROGRAM STOPS.'
write(u6,*)
write(u6,*) ' Errors occured in RASSI/RDJOB.'
call ABEND()

end subroutine RDJOB
