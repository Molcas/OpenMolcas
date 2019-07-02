************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE RDJOB(JOB,READ_STATES)
#ifdef _DMRG_
      use qcmaquis_interface_cfg
      use qcmaquis_info
#endif
      use mspt2_eigenvectors
      IMPLICIT NONE
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='RDJOB')
#include "rasdim.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "symmul.fh"
#include "rasdef.fh"
#include "rassi.fh"
#include "jobin.fh"
#include "WrkSpc.fh"
#include "Struct.fh"
#include "SysDef.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#  include "mh5.fh"
      integer :: refwfn_id

      integer :: ref_nSym, ref_lSym, ref_nBas(mxSym), ref_iSpin
      integer :: ref_nfro(mxSym), ref_nish(mxSym), ref_nrs1(mxSym),
     &           ref_nrs2(mxSym), ref_nrs3(mxSym), ref_nssh(mxSym),
     &           ref_ndel(mxSym), ref_nash(mxSym)
      integer :: ref_nactel, ref_nhole1, ref_nelec3, ref_nconf
      integer :: ref_nstates, ref_nroots
      integer, allocatable :: ref_rootid(:)

      character(1), allocatable :: typestring(:)

      real*8, allocatable :: ref_Heff(:,:), ref_energies(:)
#endif

      Real*8 Weight(MxRoot), ENUCDUMMY, AEMAX, E, HIJ
      Integer IAD, IAD15, IDISK, IERR
      Integer IPT2
      Integer ISY, IT
      Integer I, J, ISTATE, JSTATE, ISNUM, JSNUM, iAdr
      Integer LEJOB, LHEFF, NEJOB, NHEFF, NIS, NIS1, NTIT1, NMAYBE
      INTEGER JOB, NROOT0
      LOGICAL READ_STATES
#ifdef _HDF5_
      character(len=16) :: molcas_module
      character(len=8)  :: heff_string
      character(len=21) :: pt2_e_string
#endif


      CALL QENTER(ROUTINE)

#ifdef _HDF5_
************************************************************************
*
* For HDF5 formatted job files
*
************************************************************************
      If (mh5_is_hdf5(jbname(job))) Then

      IF (IPGLOB.GE.USUAL) THEN
        IF (JOB.EQ.1) THEN
          WRITE(6,*)
          WRITE(6,'(6X,80A1)') ('*',i=1,80)
          WRITE(6,'(6X,A1,78X,A1)') '*','*'
          WRITE(6,'(6X,A1,24X,A,24X,A1)')
     &       '*','     General data section     ','*'
          WRITE(6,'(6X,A1,78X,A1)') '*','*'
          WRITE(6,'(6X,80A1)') ('*',i=1,80)
        END IF
        WRITE(6,*)
        WRITE(6,*)'  Specific data for HDF5 file ',JBNAME(JOB)
        WRITE(6,*)'  -------------------------------------'
      END IF

      refwfn_id = mh5_open_file_r(jbname(job))

      call mh5_fetch_attr (refwfn_id,'MOLCAS_MODULE', molcas_module)
      call mh5_fetch_attr (refwfn_id,'SPINMULT', ref_iSpin)
      call mh5_fetch_attr (refwfn_id,'NSYM', ref_nSym)
      call mh5_fetch_attr (refwfn_id,'LSYM', ref_lSym)
      call mh5_fetch_attr (refwfn_id,'NBAS', ref_nBas)

      call mh5_fetch_attr (refwfn_id,'NACTEL', ref_nactel)
      call mh5_fetch_attr (refwfn_id,'NHOLE1', ref_nhole1)
      call mh5_fetch_attr (refwfn_id,'NELEC3', ref_nelec3)
      call mh5_fetch_attr (refwfn_id,'NCONF',  ref_nconf)
      call mh5_fetch_attr (refwfn_id,'NSTATES', ref_nstates)
      If (mh5_exists_dset(refwfn_id, 'NROOTS')) Then
        call mh5_fetch_attr (refwfn_id,'NROOTS', ref_nroots)
      Else
        ref_nroots = ref_nstates
      End If

      call mma_allocate (typestring, sum(ref_nbas(1:ref_nsym)))
      call mh5_fetch_dset (refwfn_id, 'MO_TYPEINDICES', typestring)
      call tpstr2orb (ref_nsym,ref_nbas,typestring,
     &                ref_nfro,ref_nish,ref_nrs1,ref_nrs2,ref_nrs3,
     &                ref_nssh,ref_ndel)
      ref_nash = ref_nrs1 + ref_nrs2 + ref_nrs3
      call mma_deallocate (typestring)

#ifdef _DMRG_
      If (.not.mh5_exists_dset(refwfn_id, 'CI_VECTORS').and.
     &    .not.doDMRG) Then
* Leon: TODO: This must be also extended for other DMRG interfaces
* than QCMaquis
#else
      If (.not.mh5_exists_dset(refwfn_id, 'CI_VECTORS')) then
#endif
        Write(6,'(1X,A)') 'The HDF5 file does not contain CI vectors,'
        Write(6,'(1X,A)') 'make sure it was created by rasscf/caspt2.'
        Call AbEnd()
      End If
      If (.not.mh5_exists_dset(refwfn_id, 'MO_VECTORS')) Then
        Write(6,'(1X,A)') 'The HDF5 file does not contain MO vectors,'
        Write(6,'(1X,A)') 'make sure it was created by '//
     &                    'rasscf/caspt2/nevpt2.'
        Call AbEnd()
      End If

      call mh5_fetch_attr (refwfn_id,'L2ACT', L2ACT)
      call mh5_fetch_attr (refwfn_id,'A2LEV', LEVEL)

      call mma_allocate(ref_rootid,ref_nstates)
      call mh5_fetch_attr (refwfn_id,'STATE_ROOTID', ref_rootid)
      if (read_states) then
*  Do not update the state number here, because it's already read in
*  rdjob_nstates()
*        NSTAT(JOB)=ref_nstates
*        NSTATE=NSTATE+ref_nstates
* store the root IDs of each state
        DO I=0,NSTAT(JOB)-1
          iWork(lLROOT+ISTAT(JOB)-1+I)=ref_rootid(I+1)
          iWork(lJBNUM+ISTAT(JOB)-1+I)=JOB
        END DO
      end if
      LROT1=ref_nroots
      DO I=0,NSTAT(JOB)-1
        NROOT0=iWork(lLROOT+ISTAT(JOB)-1+I)
        IF (NROOT0.GT.LROT1) THEN
          GOTO 9002
        END IF
      END DO

      if(qdpt2sc.and.(trim(molcas_module(1:6)).eq.'NEVPT2'))then
        heff_string     = 'H_EFF_SC'
        pt2_e_string    = 'STATE_PT2_ENERGIES_SC'
      else
        heff_string     = 'H_EFF'
        pt2_e_string    = 'STATE_PT2_ENERGIES'
      end if

* read the ms-caspt2/qd-nevpt2 effective hamiltonian if it is available
      If (.not.ifejob.and.mh5_exists_dset(refwfn_id, heff_string)) Then
        HAVE_HEFF=.TRUE.
        call mma_allocate(ref_Heff,ref_nstates,ref_nstates)
        call mh5_fetch_dset_array_real(refwfn_id,heff_string,ref_Heff)
        write(6,'(2x,a)')
     & ' Effective Hamiltonian from MRPT2 in action'
        write(6,'(2x,a)')
     & ' ------------------------------------------'
        DO I=1,NSTAT(JOB)
          ISTATE=ISTAT(JOB)-1+I
          DO J=1,NSTAT(JOB)
            JSTATE=ISTAT(JOB)-1+J
            iadr=(istate-1)*nstate+jstate-1
            Work(l_heff+iadr)=ref_Heff(I,J)
!           write(6,*) 'readin: Heff(',istate,',',jstate,') = ',
!    &      Work(l_heff+iadr)
!           call xflush(6)
          END DO
        END DO
        call mma_deallocate(ref_Heff)
* read the caspt2/qdnevpt2 reference energies if available
      Else If (mh5_exists_dset(refwfn_id, pt2_e_string)) Then
        HAVE_DIAG=.TRUE.
        call mma_allocate(ref_energies,ref_nstates)
        call mh5_fetch_dset_array_real(refwfn_id,
     &         pt2_e_string,ref_energies)
        DO I=1,NSTAT(JOB)
          ISTATE=ISTAT(JOB)-1+I
          Work(LREFENE+istate-1)=ref_energies(I)
        END DO
        call mma_deallocate(ref_energies)
* read rasscf energies
      Else If (mh5_exists_dset(refwfn_id, 'ROOT_ENERGIES')) Then
        HAVE_DIAG=.TRUE.
        call mma_allocate(ref_energies,ref_nroots)
        call mh5_fetch_dset_array_real(refwfn_id,
     &         'ROOT_ENERGIES',ref_energies)
        DO I=1,NSTAT(JOB)
          ISTATE=ISTAT(JOB)-1+I
          Work(LREFENE+istate-1)=ref_energies(iWork(lLROOT+ISTATE-1))
        END DO
        call mma_deallocate(ref_energies)
      End If

!     write(6,*) 'job --> ',job, 'doDMRG and doMPSSICheckpoints ',
!    & doDMRG,doMPSSICheckpoints
#ifdef _DMRG_
      ! Leon 5/12/2016: Fetch QCMaquis checkpoint names if requested
      if (doDMRG.and.doMPSSICheckpoints) then
        if(mh5_exists_dset(refwfn_id, 'QCMAQUIS_CHECKPOINT')) then
!         Write(6,'(A)') 'Reading QCMaquis checkpoint names '//
!    &    'from HDF5 files'
!         Write(6,'(A)') 'State    Checkpoint name'

          !> allocate space for the file name strings of job JOB
          call qcmaquis_info_init(job,nstat(job),1)

          DO I=1,NSTAT(JOB)
            ISTATE=ISTAT(JOB)-1+I
            call mh5_fetch_dset_array_str(refwfn_id,
     &                                    'QCMAQUIS_CHECKPOINT',
     &                                     qcm_group_names(job)
     &                                     %states(i),
     &                                     [1],
     &                                     [iWork(lLROOT+ISTATE-1)-1]
     &                                    )
!           Write(6,'(I3,A,A)') ISTATE, '   ',
!    &      trim(qcm_group_names(job)%states(i))
          END DO
        else
          call WarningMessage(2,'QCMaquis checkpoint names not found'//
     &    ' on HDF5 files. Make sure you created them with the'//
     &    ' MOLCAS version which supports them')
          call Quit_OnUserError
        end if
      end if
#endif
      if (ref_nsym.ne.nsym) then
        call WarningMessage(2,'NSYM not consistent with RunFile')
        call Quit_OnUserError
      end if
      do i=1,nsym
        if (ref_nbas(i).ne.nbasf(i)) then
          call WarningMessage(2,'NBAS not consistent with RunFile')
          call Quit_OnUserError
        end if
      end do

      NACTE(JOB)=ref_nactel
      NHOLE1(JOB)=ref_nhole1
      NELE3(JOB)=ref_nelec3
      MLTPLT(JOB)=ref_iSpin
      IRREP(JOB)=ref_lSym
      NCONF(JOB)=ref_nConf
      NROOTS(JOB)=ref_nroots

      if (job.eq.1) then
* first wavefunction file, set global variables
        DO I=1,NSYM
          NFRO(I)=0
          NISH(I)=ref_nfro(I)+ref_nish(I)
          NASH(I)=ref_nash(I)
          NRS1(I)=ref_NRS1(I)
          NRS2(I)=ref_NRS2(I)
          NRS3(I)=ref_NRS3(I)
          NOSH(I)=NISH(I)+NASH(I)
          NDEL(I)=0
          NSSH(I)=NBASF(I)-NFRO(I)-NISH(I)-NASH(I)-NDEL(I)
        END DO
      else
* subsequent wavefunction file, check against global variables
        if ( ref_nhole1.ne.nhole1(1) .or.
     &       ref_nelec3.ne.nele3(1)) then
          call WarningMessage(2,'inconsistent RAS holes/electrons')
          call Quit_OnUserError
        end if
        do i=1,nsym
          if ((ref_nfro(i)+ref_nish(i).ne.nish(i)) .or.
     &        (ref_nash(i).ne.nash(i)) .or.
     &        (ref_nrs1(i).ne.nrs1(i)) .or.
     &        (ref_nrs2(i).ne.nrs2(i)) .or.
     &        (ref_nrs3(i).ne.nrs3(i)) ) then
            call WarningMessage(2,'inconsistent orbital partitioning')
            call Quit_OnUserError
          end if
        end do
      end if

      WFTYPE='GENERAL '
      IF(ref_nactel.EQ.2*SUM(NASH(1:NSYM))) WFTYPE='CLOSED  '
      IF(ref_nactel.EQ.0) WFTYPE='EMPTY   '
      RASTYP(JOB)=WFTYPE

      IF (IPGLOB.GE.USUAL) THEN
        WRITE(6,*)'  STATE IRREP:        ',IRREP(JOB)
        WRITE(6,*)'  SPIN MULTIPLICITY:  ',MLTPLT(JOB)
        WRITE(6,*)'  ACTIVE ELECTRONS:   ',NACTE(JOB)
        WRITE(6,*)'  MAX RAS1 HOLES:     ',NHOLE1(JOB)
        WRITE(6,*)'  MAX RAS3 ELECTRONS: ',NELE3(JOB)
        WRITE(6,*)'  NR OF CONFIG:       ',NCONF(JOB)
      END IF
      IF(IPGLOB.GE.VERBOSE)
     &          WRITE(6,*)'  Wave function type WFTYPE=',WFTYPE

      call mma_deallocate(ref_rootid)
      call mh5_close_file(refwfn_id)

      Else
#endif
************************************************************************
*
* For JOBIPH/JOBMIX formatted job files
*
************************************************************************
#ifdef _DMRG_
      if (doDMRG) then
        if (doMPSSICheckpoints) then
          call WarningMessage(3, "QCMaquis checkpoint names from "//
     &   "JobIph requested. This works only with HDF5 JobIph files."//
     &   " Please make sure you use a .h5 file as JOBxxx.")
          call abend()
        else
          call WarningMessage(2, "Using old-style JobIph with DMRG "//
     &      "and hence default naming convention for checkpoint files")
        end if
      end if
#endif
      IF (IPGLOB.GE.USUAL) THEN
        IF (JOB.EQ.1) THEN
          WRITE(6,*)
          WRITE(6,'(6X,80A1)') ('*',i=1,80)
          WRITE(6,'(6X,A1,78X,A1)') '*','*'
          WRITE(6,'(6X,A1,24X,A,24X,A1)')
     &       '*','     General data section     ','*'
          WRITE(6,'(6X,A1,78X,A1)') '*','*'
          WRITE(6,'(6X,80A1)') ('*',i=1,80)
        END IF
        WRITE(6,*)
        WRITE(6,*)'  Specific data for JOBIPH file ',JBNAME(JOB)
        WRITE(6,*)'  -------------------------------------'
      END IF
C Open JOBIPH file:
      CALL DANAME(LUIPH,JBNAME(JOB))
C READ TABLE OF CONTENTS ON THIS JOBIPH FILE:
      IAD=0
      CALL IDAFILE(LUIPH,2,ITOC15,30,IAD)
C SCATTER-READ VARIOUS DATA:
* PAM Mar2014: Note that ENUC1 (=POTNUC) replaced by dummy placeholder
      IAD=ITOC15(1)
      Call WR_RASSCF_Info(LUIPH,2,IAD,
     &                    NACTE1,MPLET1,NSYM1,LSYM1,
     &                    NFRO1,NISH1,NASH1,NDEL1,NBAS1,mxSym,
     &                    NAME,LENIN8*mxOrb,NCONF1,HEAD1,2*72,
     &                    TITLE1,4*mxTit*18,
     &                    ENUCDUMMY,LROT1,NROOT1,
     &                    IROOT1,mxRoot,NRS11,NRS21,NRS31,
     &                    NHOL11,NELE31,IPT2,Weight)
C Response field contribution to zero-electron energies
C is added in GETH1.
      IF (READ_STATES) THEN
* Do not update the state number here, because it's already read in
* rdjob_nstates()
!        ISTAT(JOB)=NSTATE+1
!        NSTAT(JOB)=NROOT1
!        NSTATE=NSTATE+NROOT1
* store the root IDs of each state

*If unset yet, set now
        If (iWork(lLROOT+ISTAT(JOB)-1).eq.0) Then
          DO I=0,NSTAT(JOB)-1
            iWork(lLROOT+ISTAT(JOB)-1+I)=IROOT1(I+1)
            iWork(lJBNUM+ISTAT(JOB)-1+I)=JOB
          End DO
        End If
      END IF
      DO I=0,NSTAT(JOB)-1
        NROOT0=iWork(lLROOT+ISTAT(JOB)-1+I)
        IF (NROOT0.GT.LROT1) THEN
          GOTO 9002
        END IF
      END DO

C Using energy data from JobIph?
      IF(IFEJOB) THEN
        NEJOB=MXROOT*MXITER
        CALL GETMEM('EJOB','ALLO','REAL',LEJOB,NEJOB)
        IAD=ITOC15(6)
        CALL DDAFILE(LUIPH,2,WORK(LEJOB),NEJOB,IAD)
C Note that there is no info on nr of iterations
C so we cannot know what energies to pick...
C Let us make a guess: The correct set of energy values in the
C table of energies/iteration is the last one with not all zeroes.
        NMAYBE=0
        DO IT=1,MXITER
          AEMAX=0.0D0
          DO I=1,MXROOT
            E=WORK(LEJOB+MXROOT*(IT-1)+(I-1))
            AEMAX=MAX(AEMAX,ABS(E))
          END DO
          IF(ABS(AEMAX).LE.1.0D-12) GOTO 11
          NMAYBE=IT
        END DO
  11    CONTINUE
        IF(NMAYBE.EQ.0) THEN
          WRITE(6,*)' Sorry. Keyword ''EJOB'' has been used'
          WRITE(6,*)' but there are no energies available on'
          WRITE(6,*)' the JOBIPH file nr', JOB
          CALL ABEND()
        END IF
        HAVE_DIAG=.TRUE.

C Put these energies into diagonal of Hamiltonian:
        DO I=1,NSTAT(JOB)
          ISTATE=ISTAT(JOB)-1+I
#ifdef _DMRG_
          if (doDMRG) then
            E=WORK(LEJOB-1+iWork(lLROOT+ISTATE-1)
     &        -ISTAT(JOB)+1+MXROOT*(NMAYBE-1))
          else
#endif
          E=WORK(LEJOB-1+iWork(lLROOT+ISTATE-1)+MXROOT*(NMAYBE-1))
#ifdef _DMRG_
          endif
#endif
          Work(LREFENE+istate-1)=E
        END DO
        CALL GETMEM('EJOB','FREE','REAL',LEJOB,NEJOB)
      END IF
C Using effective Hamiltonian from JobIph file?
      IF(IFHEFF) THEN
        IF(ITOC15(15).NE.-1) THEN
          WRITE(6,*)'RDJOB Error: HEFF not found on JOBIPH.'
          WRITE(6,*)'The HEFF keyword was used, but the JOBIPH file'
          WRITE(6,*)'uses an old layout where this data field is'
          WRITE(6,*)'not present. Recompute JOBIPH file, or put'
          WRITE(6,*)'effective Hamiltonian in input after keyword'
          WRITE(6,*)'HEXT. Program stops here.'
          CALL ABEND()
        END IF
        HAVE_HEFF=.TRUE.
        NHEFF=LROT1**2
        CALL GETMEM('HEFF','ALLO','REAL',LHEFF,NHEFF)
        IAD15=ITOC15(17)
        CALL DDAFILE(LUIPH,2,WORK(LHEFF),NHEFF,IAD15)
        DO I=1,NSTAT(JOB)
          ISTATE=ISTAT(JOB)-1+I
          ISNUM=iWork(lLROOT+ISTATE-1)
          DO J=1,NSTAT(JOB)
            JSTATE=ISTAT(JOB)-1+J
            JSNUM=iWork(lLROOT+JSTATE-1)
            HIJ=WORK(LHEFF-1+ISNUM+LROT1*(JSNUM-1))
            iadr=(istate-1)*nstate+jstate-1
            Work(l_heff+iadr)=HIJ
          END DO
        END DO
        CALL GETMEM('HEFF','FREE','REAL',LHEFF,NHEFF)
      END IF
C Read the level to orbital translations
      IDISK=ITOC15(18)
      CALL IDAFILE(LUIPH,2,L2ACT,MXLEV,IDISK)
      CALL IDAFILE(LUIPH,2,LEVEL,MXLEV,IDISK)
C Close JobIph file
      CALL DACLOS(LUIPH)

C The RASSCF program is not certain to give consistent data. For
C pure CASSCF cases, it may not bother to set the NRS1..NRS3 arrays.
C Check and repair:
      IF(NHOL11+NELE31.EQ.0) THEN
        IERR=0
        DO I=1,NSYM1
          IF(NRS11(I).NE.0) IERR=1
          IF(NRS21(I).NE.NASH1(I)) IERR=1
          IF(NRS31(I).NE.0) IERR=1
        END DO
        IF(IERR.EQ.1) THEN
          WRITE(6,*)
          WRITE(6,*)' (NOTE: The nr of RAS1, RAS2 and RAS3 orbitals'//
     &              ' as recorded on the JOBIPH file do not match'
          WRITE(6,*)' the number of active orbitals. But this is a'//
     &              ' pure CASSCF case. Maybe the RASSCF programmer did'
          WRITE(6,*)' not bother with the RAS1..RAS3 arrays in that'//
     &              ' case. RASSI will reset these arrays as needed.)'
          WRITE(6,*)
          DO I=1,NSYM1
            NRS11(I)=0
            NRS21(I)=NASH1(I)
            NRS31(I)=0
          END DO
        END IF
      END IF

      IF(JOB.EQ.1) THEN
C FIRST JOB FILE. TRANSFER DATA TO COMMON:
        NSYM=NSYM1
        DO I=1,NSYM
          NFRO(I)=0
          NISH(I)=NFRO1(I)+NISH1(I)
          NASH(I)=NASH1(I)
          NRS1(I)=NRS11(I)
          NRS2(I)=NRS21(I)
          NRS3(I)=NRS31(I)
          NOSH(I)=NISH(I)+NASH(I)
          NDEL(I)=0
          NBASF(I)=NBAS1(I)
          NSSH(I)=NBASF(I)-NFRO(I)-NISH(I)-NASH(I)-NDEL(I)
        END DO
      ELSE
C THIS IS NOT THE FIRST JOBIPH.
C CHECK THAT DATA IS CONSISTENT WITH EARLIER:
        IF(NSYM1.NE.NSYM) GOTO 9001
        IF(NHOL11.ne.NHOLE1(JOB-1)) GOTO 9003
        IF(NELE31.ne.NELE3(JOB-1))  GOTO 9003
        DO ISY=1,NSYM1
          NIS1=NISH1(ISY)+NFRO1(ISY)
          NIS =NISH (ISY)
          IF(NIS1.NE.NIS ) GOTO 9004
          IF(NRS11(ISY).NE.NRS1 (ISY)) GOTO 9005
          IF(NRS21(ISY).NE.NRS2 (ISY)) GOTO 9005
          IF(NRS31(ISY).NE.NRS3 (ISY)) GOTO 9005
          IF(NBAS1(ISY).NE.NBASF(ISY)) GOTO 9006
        END DO
      END IF

C DATA PARTICULAR TO THIS JOBIPH:
      IF (IPGLOB.GE.USUAL) THEN
        WRITE(6,*)
        WRITE(6,*)'  Header from SEWARD:'
        WRITE(6,'(7X,36A2)')(HEAD1(I),I=1,36)
        WRITE(6,'(7X,36A2)')(HEAD1(I),I=37,72)
C NOTE: AT PRESENT, JOBIPH FILE GIVES NO INFORMATION ON THE
C AMOUNT OF TITLE LINES.
        NTIT1=1
        WRITE(6,*)
        WRITE(6,*)'  CASSCF title (first line only):'
        WRITE(6,'(7X,18A4)')((TITLE1(I,J),I=1,18),J=1,NTIT1)
        WRITE(6,*)
        WRITE(6,*)'  STATE IRREP:        ',LSYM1
        WRITE(6,*)'  SPIN MULTIPLICITY:  ',MPLET1
        WRITE(6,*)'  ACTIVE ELECTRONS:   ',NACTE1
        WRITE(6,*)'  MAX RAS1 HOLES:     ',NHOL11
        WRITE(6,*)'  MAX RAS3 ELECTRONS: ',NELE31
        WRITE(6,*)'  NR OF CONFIG:       ',NCONF1
      END IF
      WFTYPE='GENERAL '
*      IF(MPLET1.EQ.(SUM(NASH(1:NSYM))+1)) WFTYPE='HISPIN  '
* Note: the HISPIN case may be buggy and is not used presently.
      IF(MPLET1.EQ.(SUM(NASH(1:NSYM))+1)) THEN
       write(6,*)' This wave function is of HISPIN type.'
       write(6,*)' However, the special handling for that case'
       write(6,*)' is suspected to be buggy. So the variable'
       write(6,*)' WFTYPE is set to GENERAL.'
      END IF
      IF(NACTE1.EQ.2*SUM(NASH(1:NSYM))) WFTYPE='CLOSED  '
      IF(NACTE1.EQ.0) WFTYPE='EMPTY   '
      RASTYP(JOB)=WFTYPE
      IF(IPGLOB.GE.VERBOSE)
     &          WRITE(6,*)'  Wave function type WFTYPE=',WFTYPE
      NACTE(JOB)=NACTE1
      MLTPLT(JOB)=MPLET1
      IRREP(JOB)=LSYM1
      NHOLE1(JOB)=NHOL11
      NELE3(JOB)=NELE31
      NCONF(JOB)=NCONF1
C Where is the CMO data set stored?
      IDCMO(JOB)=ITOC15(2)
      IF(IPT2.NE.0) IDCMO(JOB)=ITOC15(9)

#ifdef _HDF5_
      End If
#endif

      CALL XFLUSH(6)
      CALL QEXIT(ROUTINE)
      RETURN
************************************************************************
*
* Error exits
*
************************************************************************
9001  WRITE(6,*) ' SYMMETRY GROUPS MUST BE EQUAL.'
      WRITE(6,*) ' NSYM1:',NSYM1,'NSYM :',NSYM
      GOTO 9010
9002  WRITE(6,*) ' ROOT NOT AVAILABLE.'
      WRITE(6,*) '             REQUESTED ROOT:',NROOT0
      WRITE(6,*) '  MAXIMUM ROOT IN THIS FILE:',LROT1
      GOTO 9010
9003  WRITE(6,*) ' RAS SPECIFICATIONS DIFFER.'
      WRITE(6,*) '     THIS STATE: MAX NR OF RAS-1 HOLES:',NHOL11
      WRITE(6,*) '             MAX NR OF RAS-3 ELECTRONS:',NELE31
      WRITE(6,*) ' PREVIOUS STATE: MAX NR OF RAS-1 HOLES:',NHOLE1(JOB-1)
      WRITE(6,*) '             MAX NR OF RAS-3 ELECTRONS:',NELE3(JOB-1)
      GOTO 9010
9004  WRITE(6,*) ' NR. OF (FROZEN+INACTIVE) ORBITALS DIFFER.'
      WRITE(6,'(A,8I4)')' THIS STATE:',(NFRO1(I)+NISH1(I),I=1,NSYM1)
      WRITE(6,'(A,8I4)')'   PREVIOUS:',(NISH(I),I=1,NSYM )
      GOTO 9010
9005  WRITE(6,*)' NR. OF ACTIVE ORBITALS DIFFER.'
      WRITE(6,'(A,8I4)')' THIS STATE, ACTIVE:',(NASH1(I),I=1,NSYM1)
      WRITE(6,'(A,8I4)')'               RAS1:',(NRS11(I),I=1,NSYM1)
      WRITE(6,'(A,8I4)')'               RAS2:',(NRS21(I),I=1,NSYM1)
      WRITE(6,'(A,8I4)')'               RAS3:',(NRS31(I),I=1,NSYM1)
      WRITE(6,'(A,8I4)')'   PREVIOUS, ACTIVE:',(NASH (I),I=1,NSYM )
      WRITE(6,'(A,8I4)')'               RAS1:',(NRS1 (I),I=1,NSYM1)
      WRITE(6,'(A,8I4)')'               RAS2:',(NRS2 (I),I=1,NSYM1)
      WRITE(6,'(A,8I4)')'               RAS3:',(NRS3 (I),I=1,NSYM1)
      GOTO 9010
9006  WRITE(6,*)' NR. OF BASIS FUNCTION DIFFER.'
      WRITE(6,'(A,8I4)')' THIS JOBIPH:',(NBAS1(I),I=1,NSYM1)
      WRITE(6,'(A,8I4)')'    PREVIOUS:',(NBASF(I),I=1,NSYM )
9010  CONTINUE
      WRITE(6,*)' DATA IN JOBIPH FILE NAMED ',TRIM(JBNAME(JOB)),' WERE'
      WRITE(6,*)' INCONSISTENT WITH EARLIER DATA. PROGRAM STOPS.'
      WRITE(6,*)
      WRITE(6,*)' Errors occured in RASSI/RDJOB.'
      CALL XFLUSH(6)
      CALL ABEND()

      END

************************************************************************
*                                                                      *
*     Only read the number of states                                   *
*                                                                      *
************************************************************************
      Subroutine rdjob_nstates(JOB)
      IMPLICIT NONE
#include "rasdim.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "jobin.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#  include "mh5.fh"
      integer :: refwfn_id
      integer :: ref_nstates
#endif
      Real*8 Weight(MxRoot), ENUCDUMMY
      Integer job,iad,ipt2
************************************************************************
*
* For HDF5 formatted job files
*
************************************************************************
#ifdef _HDF5_
      If (mh5_is_hdf5(jbname(job))) Then
        refwfn_id = mh5_open_file_r(jbname(job))
        call mh5_fetch_attr (refwfn_id,'NSTATES', ref_nstates)
* update the state offset, number of states, and total number of states
        ISTAT(JOB)=NSTATE+1
        NSTAT(JOB)=ref_nstates
        NSTATE=NSTATE+ref_nstates
        call mh5_close_file(refwfn_id)
      Else
#endif
      CALL DANAME(LUIPH,JBNAME(JOB))
C READ TABLE OF CONTENTS ON THIS JOBIPH FILE:
      IAD=0
      CALL IDAFILE(LUIPH,2,ITOC15,30,IAD)
C SCATTER-READ VARIOUS DATA:
      IAD=ITOC15(1)
      Call WR_RASSCF_Info(LUIPH,2,IAD,
     &                    NACTE1,MPLET1,NSYM1,LSYM1,
     &                    NFRO1,NISH1,NASH1,NDEL1,NBAS1,mxSym,
     &                    NAME,LENIN8*mxOrb,NCONF1,HEAD1,2*72,
     &                    TITLE1,4*mxTit*18,
     &                    ENUCDUMMY,LROT1,NROOT1,
     &                    IROOT1,mxRoot,NRS11,NRS21,NRS31,
     &                    NHOL11,NELE31,IPT2,Weight)
* update the state offset, number of states, and total number of states
      ISTAT(JOB)=NSTATE+1
      NSTAT(JOB)=NROOT1
      NSTATE=NSTATE+NROOT1
      CALL DACLOS(LUIPH)
#ifdef _HDF5_
      EndIf
#endif
      end
