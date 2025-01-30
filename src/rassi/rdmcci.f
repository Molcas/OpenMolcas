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
      use rassi_aux, only: ipglob
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NJOB, NCONF1, MINAME
      use cntrl, only: LuMck

      IMPLICIT None
C Purpose: Read in the derivatives of CI array derivatives
C from MCKINT file, with respect to some displacement IDISP.
C ISYMP is the symmetry irrep label of the derivatives.
      Integer JOB, IDISP, ISYMP, nArray
      CHARACTER(LEN=8) LABEL
      Real*8 ARRAY(NARRAY)

      Real*8, Allocatable:: TEMP(:)
      Integer IRC, IOPT, NTEMP, ISCODE

      IF(JOB.LT.1 .OR. JOB.GT.NJOB) THEN
        WRITE(6,*)' RDMCI: Invalid JOB parameter.'
        WRITE(6,*)' JOB:',JOB
        CALL ABEND()
      END IF

      IF(IPGLOB.GE.3) THEN
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
      CALL mma_allocate(TEMP,NTEMP,Label='TEMP')
C Read 1-electron integral derivatives:
      IRC=NTEMP
      CALL dRDMCK(IRC,IOPT,LABEL,IDISP,TEMP,ISCODE)
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
      CALL DCOPY_(NTEMP,TEMP,1,ARRAY,1)
C Get rid of temporary buffer
      CALL mma_deallocate(TEMP)

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

      END SUBROUTINE RDMCCI
