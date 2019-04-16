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
      SUBROUTINE RDMCCI(JOB,IDISP,LABEL,ISYMP,NARRAY,ARRAY)
      IMPLICIT REAL*8 (A-H,O-Z)
C Purpose: Read in the derivatives of CI array derivatives
C from MCKINT file, with respect to some displacement IDISP.
C ISYMP is the symmetry irrep label of the derivatives.
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='RDMCCI')
#include "rasdim.fh"
#include "SysDef.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "rassi.fh"
#include "jobin.fh"
#include "Struct.fh"
#include "symmul.fh"
#include "WrkSpc.fh"
      DIMENSION ARRAY(NARRAY)
      CHARACTER*8 LABEL

      CALL QENTER(ROUTINE)

      IF(JOB.LT.1 .OR. JOB.GT.NJOB) THEN
        WRITE(6,*)' RDMCI: Invalid JOB parameter.'
        WRITE(6,*)' JOB:',JOB
        CALL ABEND()
      END IF

      IF(IPGLOB.GE.VERBOSE) THEN
        WRITE(6,*)' RDMCCI called for JOB=',JOB
        WRITE(6,*)' perturbed by displacement nr.',IDISP
        WRITE(6,*)' MckInt file name:',MINAME(JOB)
        WRITE(6,*)' Irrep label   ISYMP=',ISYMP
        WRITE(6,*)' Length NARRAY=',NARRAY
      END IF

C Open MCKINT file:
      IRC=-1
      IOPT=0
      CALL OPNMCK(IRC,IOPT,MINAME(JOB),LUMCK)
      IF(IRC.NE.0) THEN
        WRITE(6,*)'RDMCCI: Failed to open '//MINAME(JOB)
        WRITE(6,*)'Unit nr LUMCK=',LUMCK
        WRITE(6,*)'Option code IOPT=',IOPT
        WRITE(6,*)'Return code IRC =',IRC
        CALL ABEND()
      END IF

C Read MCKINT file:
      NTEMP=NCONF1
      IOPT=0
      ISCODE=2**(ISYMP-1)
C Get temporary buffer to read data by RDMCK calls
      CALL GETMEM('RDMCCI','ALLO','REAL',LTEMP,NTEMP)
C Read 1-electron integral derivatives:
      IRC=NTEMP
      CALL dRDMCK(IRC,IOPT,LABEL,IDISP,WORK(LTEMP),ISCODE)
      IF(IRC.NE.0) THEN
        WRITE(6,*)'RDMCCI: RDMCCI failed to read '//MINAME(JOB)
        WRITE(6,*)'  Displacement IDISP=',IDISP
        WRITE(6,*)'    Option code IOPT=',IOPT
        WRITE(6,*)'Symmetry code ISCODE=',ISCODE
        WRITE(6,*)'    Return code IRC =',IRC
        CALL ABEND()
      END IF

      IF(NTEMP.GT.NARRAY) THEN
        WRITE(6,*)'RDMCCI: Output ARRAY has insufficient length.'
        WRITE(6,*)' Input parameter NARRAY=',NARRAY
        WRITE(6,*)' Needed size       NTEMP=',NTEMP
        CALL ABEND()
      END IF
C Move buffer integrals into ARRAY in proper format:
      CALL DCOPY_(NTEMP,WORK(LTEMP),1,ARRAY,1)
C Get rid of temporary buffer
      CALL GETMEM('RDMCCI','FREE','REAL',LTEMP,NTEMP)

C Close MCKINT file:
      IRC=-1
      IOPT=0
      CALL CLSMCK(IRC,IOPT)
      IF(IRC.NE.0) THEN
        WRITE(6,*)'RDMCCI: Failed to close '//MINAME(JOB)
        WRITE(6,*)'Unit nr LUMCK=',LUMCK
        WRITE(6,*)'Option code IOPT=',IOPT
        WRITE(6,*)'Return code IRC =',IRC
        CALL ABEND()
      END IF

      CALL QEXIT(ROUTINE)
      RETURN
      END
