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
      SUBROUTINE INPCTL_RASSI()
      use rassi_global_arrays, only: HAM, ESHFT, HDIAG, JBNUM, LROOT
#ifdef _DMRG_
      use rasscf_global, only: doDMRG
      use qcmaquis_interface_cfg
      use qcmaquis_info, only: qcmaquis_info_init, qcm_prefixes
      use qcmaquis_interface_mpssi, only: qcmaquis_mpssi_init
      use cntrl, only: NACTE
      use rassi_data, only: NASH
#endif
      use mspt2_eigenvectors
      use stdalloc, only: mma_allocate, mma_deallocate
      use cntrl, only: RefEne, HEff
      use Cntrl, only:  NSTATE, NJOB, IFHEXT, IFShft, IfHDia, ISTAT,
     &                  MLTPLT, NSTAT, MXJOB
      use cntrl, only: ATLBL, IGROUP, nAtoms, nGroup
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: ENUC,NBASF
      IMPLICIT NONE

      LOGICAL READ_STATES
      INTEGER JOB, i


* get basic info from runfile
      Call Get_iArray('nBas',nBasF,nSym)
      Call Get_dscalar('PotNuc',ENUC)

C Read data from the ONEINT file:
      CALL GETCNT(NGROUP,IGROUP,NATOMS,ATLBL)

      NSTATE=0
C Read (and do some checking) the standard input.
      CALL READIN_RASSI
* if there have been no states selected at this point, we need to read
* the states later from the job files.
      IF(NSTATE.EQ.0) THEN
        READ_STATES=.TRUE.
        DO JOB=1,NJOB
          call rdjob_nstates(JOB)
        END DO
* store the root IDs of each state
        Call mma_allocate(JBNUM,nState,Label='JBNUM')
        Call mma_allocate(LROOT,nState,Label='LROOT')
        LROOT(:)=0
        Do JOB=1,NJOB
          DO I=0,NSTAT(JOB)-1
            JBNUM(ISTAT(JOB)+I)=JOB
          End Do
        End Do
      ELSE
        READ_STATES=.FALSE.
      END IF

#ifdef _DMRG_
      !> initialize DMRG interface
      if (doDMRG) then
        !> initialize only the qcm file name array (one for each job) and initialize the DMRG interface later
        call qcmaquis_info_init(njob,-1,0)
      endif
#endif

      !> initialize eigenvector array for mspt2 hamiltonians
      call init_mspt2_eigenvectors(njob,-1,0)
* Allocate a bunch of stuff
      Call mma_allocate(REFENE,NSTATE,Label='RefEne')
      Call mma_allocate(HEFF,NSTATE,NSTATE,Label='HEff')
      HEff(:,:)=0.0D0
      If (.not.IFHEXT) Then
        Call mma_allocate(HAM,nState,nState,Label='HAM')
        HAM(:,:)=0.0D0
      EndIf
      If (.not.IFSHFT) Then
        Call mma_allocate(ESHFT,nSTATE,Label='ESHFT')
        ESHFT(:)=0.0D0
      EndIf
      If (.not.IFHDIA) Call mma_Allocate(HDIAG,nState,Label='HDIAG')

C Read information on the job files and check for consistency
      DO JOB=1,NJOB
        CALL RDJOB(JOB,READ_STATES)
      END DO

C Number of active orbitals is taken from the first JobIph. MPS-SI cannot
C handle different active spaces per JobIph, but this is checked elsewhere
#ifdef _DMRG_
      if (doDMRG)then
        !> stupid info.h defines "sum", so I cannot use the intrinsic sum function here...
        qcmaquis_param%L = 0; do i = 1, nsym; qcmaquis_param%L =
     &  qcmaquis_param%L + nash(i); end do

        ! Initialise the new MPSSI interface
        call qcmaquis_mpssi_init(qcm_prefixes,
     &                           LROOT,NSTAT(1),NJOB)

      ! Check if number of active electrons is the same for all job files
      ! Otherwise, quit on error, as Dyson orbitals are not supported yet with DMRG
      if (NJOB.gt.1) then
        JOB=NACTE(1)
        do i=2,NJOB
          if (NACTE(i).ne.JOB) then
            Call WarningMessage(2,'Number of active electrons')
            Write(6,*)' is not the same in different JOBIPH files'
            Write(6,*)' Dyson orbitals are not yet supported in MPSSI.'
            Call Quit_OnUserError()
          end if
        end do
      end if

      end if
#endif

* set orbital partitioning data
      CALL WFNSIZES_RASSI()

C Added by Ungur Liviu on 04.11.2009
C Addition of NJOB,MSJOB and MLTPLT on RunFile.

      CALL Put_iscalar('NJOB_SINGLE',NJOB)
      CALL Put_iscalar('MXJOB_SINGLE',MXJOB)
      CALL Put_iArray('MLTP_SINGLE',MLTPLT,MXJOB)

      CALL Put_iArray('NSTAT_SINGLE',NSTAT,MXJOB)
!     CALL Put_iArray('ISTAT_SINGLE',ISTAT,MXJOB)
C
C .. and print it out

C Additional input processing. Start writing report.
      CALL INPPRC()
*
      Call mma_deallocate(REFENE)
      Call mma_deallocate(HEff)
C
      END
