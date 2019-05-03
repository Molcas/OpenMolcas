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
#ifdef _DMRG_
      use qcmaquis_interface_cfg
      use qcmaquis_interface_environment,
     & only: initialize_dmrg_rassi
      use qcmaquis_info
#endif
      use mspt2_eigenvectors
      IMPLICIT NONE
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='INPCTL')
#include "WrkSpc.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "itmax.fh"
#include "info.fh"
#include "centra.fh"
#include "rasdef.fh"
#include "cntrl.fh"
#ifdef _HDF5_
#  include "mh5.fh"
#endif

      LOGICAL READ_STATES
      INTEGER JOB, i

      CALL QENTER(ROUTINE)

* get basic info from runfile
      Call Get_iScalar('nSym',nSym)
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
        Call GetMem('JBNUM','Allo','Inte',LJBNUM,NSTATE)
        Call GetMem('LROOT','Allo','Inte',LLROOT,NSTATE)
        call izero(iWork(LLROOT),NSTATE)
        Do JOB=1,NJOB
          DO I=0,NSTAT(JOB)-1
            iWork(lJBNUM+ISTAT(JOB)-1+I)=JOB
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
      Call GetMem('REFENE','Allo','Real',LREFENE,NSTATE)
      Call GetMem('HEFF','Allo','Real',L_HEFF,NSTATE**2)
      Call dzero(Work(L_HEFF),NSTATE**2)
      If (.not.IFHEXT) Then
        Call GetMem('HAM','Allo','Real',LHAM,NSTATE**2)
        call dzero(Work(LHAM),NSTATE**2)
      EndIf
      If (.not.IFSHFT) Then
        Call GetMem('ESHFT','Allo','Real',LESHFT,NSTATE)
        call dzero(Work(LESHFT),NSTATE)
      EndIf
      If (.not.IFHDIA) Call GetMem('HDIAG','Allo','Real',LHDIAG,NSTATE)

C Read information on the job files and check for consistency
      DO JOB=1,NJOB
        CALL RDJOB(JOB,READ_STATES)
      END DO

C Number of active orbitals is taken from the first JobIph. MPS-SI cannot
C handle different active spaces per JobIph, but this is checked elsewhere
#ifdef _DMRG_
      if (doDMRG)then
        !> stupid info.h defines "sum", so I cannot use the intrinsic sum function here...
        dmrg_external%norb = 0; do i = 1, nsym; dmrg_external%norb =
     &  dmrg_external%norb + nash(i); end do
        !> initialize the MPS-SI interface
        call initialize_dmrg_rassi(nstate)
      end if
#endif

* set orbital partitioning data
      CALL WFNSIZES_RASSI

C Added by Ungur Liviu on 04.11.2009
C Addition of NJOB,MSJOB and MLTPLT on RunFile.

      CALL Put_iscalar('NJOB_SINGLE',NJOB)
      CALL Put_iscalar('MXJOB_SINGLE',MXJOB)
      CALL Put_iArray('MLTP_SINGLE',MLTPLT,MXJOB)

C
C .. and print it out
CTEST      CALL PRINF()
C Set up tables of coordinates and differentiated nuclei:
#if 0
      IF(NONA) THEN
        CALL MKDISP
      END IF
#endif


C Additional input processing. Start writing report.
      CALL INPPRC
*
      Call GetMem('REFENE','Free','Real',LREFENE,NSTATE)
      Call GetMem('HEFF','Free','Real',L_HEFF,NSTATE**2)
C
      CALL QEXIT(ROUTINE)
      RETURN
      END
