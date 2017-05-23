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
      SUBROUTINE LUCIA()
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Parameters for dimensioning
#include "mxpdim.fh"
*.Memory
#include "WrkSpc.fh"
*.File numbers
#include "clunit.fh"
*.Print flags
#include "cprnt.fh"
#include "lucinp.fh"
#include "cstate.fh"
#include "crun.fh"
#include "cicisp.fh"
#include "oper.fh"
#include "cgas.fh"

#include "glbbas.fh"
#include "rasscf_lucia.fh"
#include "warnings.fh"
*.Scratch : A character line
*
      CALL QENTER('REST ')
*.    No floating point underflow
      CALL XUFLOW
*. Assign diskunits
c      IF (ENVIRO(1:6) .EQ. 'RASSCF') THEN
         CALL DISKUN2
c      ELSE
c         CALL DISKUN
c      ENDIF
*. Header
*. From shells to orbitals
      CALL ORBINF(IPRORB)
*. Number of string types
      CALL STRTYP_GAS(IPRSTR)
*. Divide orbital spaces into inactive/active/secondary
      CALL GASSPC
*. Symmetry information
      CALL SYMINF_LUCIA(IPRORB)
*. Number of integrals
      CALL INTDIM(IPRORB)
*. Static memory, initialization and  allocation
c      IF (ENVIRO(1:6) .NE. 'RASSCF') then
c         KBASE = 1
c         KADD = MXPWRD
c         CALL MEMMAN(KBASE,KADD,'INI   ',IDUMMY,'DUMMY')
c      ENDIF
      CALL ALLOC_LUCIA
*. Read in integrals
         IF(NOINT.EQ.0) THEN
            CALL INTIM
         ELSE
            WRITE(6,*) ' No integrals imported '
         END IF
*. READ in MO-AO matrix
c      IF(NOMOFL.EQ.0) CALL GET_CMOAO(WORK(KMOAOIN))
*. Internal string information (stored in WORK, bases in /STRBAS/)
c       CALL STRINF_GAS(WORK,IPRSTR)
      CALL STRINF_GAS(iWORK,IPRSTR) ! really iWORK??? in master?
*. Internal subspaces
      CALL LCISPC(IPRCIX)
*
*. Symmetry of reference
c      IF(PNTGRP.GT.1) CALL MLSM(IREFSM,IREFPA,IREFSM,'CI',1)
      IF(NOINT.EQ.1) THEN
        WRITE(6,*) ' End of calculation without integrals'
        CALL QSTAT(' ')
        CALL QUIT(_RC_ALL_IS_WELL_)
      END IF
*
*. Last space where CI vectors were stored
*
      ISTOSPC = 0
      IF(IRESTR.EQ.1) ISTOSPC = 1
      IRESTR_ORIG=IRESTR
*
      LBLOCK = MXSOOB
      LBLOCK = MAX(LBLOCK,LCSBLK)
* JESPER : Should reduce I/O
c         IF (ENVIRO(1:6).EQ.'RASSCF') THEN
      LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
      IF(PSSIGN.NE.0.0D0) LBLOCK = INT(2.0D0*XISPSM(IREFSM,1))

      CALL GETMEM('VEC1  ','ALLO','REAL',KCI_POINTER,LBLOCK)
      CALL GETMEM('VEC2  ','ALLO','REAL',KSIGMA_POINTER,LBLOCK)
      CALL QEXIT('REST ')
      RETURN
      END
