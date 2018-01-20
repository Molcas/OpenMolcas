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

      INTEGER JOB,i

      CALL QENTER(ROUTINE)

* get basic info from runfile
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBasF,nSym)
      Call Get_dscalar('PotNuc',ENUC)

C Read data from the ONEINT file:
      CALL GETCNT(NGROUP,IGROUP,NATOMS,ATLBL,COOR)

      NSTATE=0
C Read (and do some checking) the standard input.
      CALL READIN_RASSI
* if there have been no states selected at this point, we need to read
* the states from the job files.
      IF(NSTATE.EQ.0) THEN
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
      END IF

* Allocate a bunch of stuff
      Call GetMem('REFENE','Allo','Real',LREFENE,NSTATE)
      L_HEFF=ip_Dummy
      If (ifheff) Then
         Call GetMem('HEFF','Allo','Real',L_HEFF,NSTATE**2)
         Call dzero(Work(L_HEFF),NSTATE**2)
      EndIf
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
        CALL RDJOB(JOB)
      END DO
* set orbital partitioning data
      CALL WFNSIZES

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
      If (ifheff) Call GetMem('HEFF','Free','Real',L_HEFF,NSTATE**2)
C
      CALL QEXIT(ROUTINE)
      RETURN
      END
