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
      SUBROUTINE LUCIA()
      use stdalloc, only: mma_allocate
      use GLBBAS, only: CI_VEC, SIGMA_VEC
      use lucia_data, only: MXSOOB,XISPSM
      use lucia_data, only: IPRCIX,IPRORB,IPRSTR
      use lucia_data, only: NOINT,LCSBLK
      use lucia_data, only: IREFSM,PSSIGN
!
      IMPLICIT NONE
!. Parameters for dimensioning
#include "warnings.h"
!.Scratch : A character line
      Integer LBLOCK
!
!.    No floating point underflow
      !CALL XUFLOW
!. Assign diskunits
!      IF (ENVIRO(1:6) .EQ. 'RASSCF') THEN
         CALL DISKUN2()
!      ELSE
!         CALL DISKUN
!      ENDIF
!. Header
!. From shells to orbitals
      CALL ORBINF(IPRORB)
!. Number of string types
      CALL STRTYP_GAS(IPRSTR)
!. Divide orbital spaces into inactive/active/secondary
      CALL GASSPC()
!. Symmetry information
      CALL SYMINF_LUCIA(IPRORB)
!. Number of integrals
      CALL INTDIM(IPRORB)
!. Static memory, initialization and  allocation
!      IF (ENVIRO(1:6) .NE. 'RASSCF') then
!         KBASE = 1
!         KADD = MXPWRD
!         CALL MEMMAN(KBASE,KADD,'INI   ',IDUMMY,'DUMMY')
!      ENDIF
      CALL ALLOC_LUCIA()
!. Read in integrals
         IF(NOINT.EQ.0) THEN
            CALL INTIM()
         ELSE
            WRITE(6,*) ' No integrals imported '
         END IF
!. READ in MO-AO matrix
!      IF(NOMOFL.EQ.0) CALL GET_CMOAO(MOAOIN)
!. Internal string information
      CALL STRINF_GAS(IPRSTR)
!. Internal subspaces
      CALL LCISPC(IPRCIX)
!
!. Symmetry of reference
!      IF(PNTGRP.GT.1) CALL MLSM(IREFSM,IREFPA,IREFSM,'CI',1)
      IF(NOINT.EQ.1) THEN
        WRITE(6,*) ' End of calculation without integrals'
        CALL QUIT(_RC_ALL_IS_WELL_)
      END IF
!
      LBLOCK = MXSOOB
      LBLOCK = MAX(LBLOCK,LCSBLK)
! JESPER : Should reduce I/O
!         IF (ENVIRO(1:6).EQ.'RASSCF') THEN
      LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
      IF(PSSIGN.NE.0.0D0) LBLOCK = INT(2.0D0*XISPSM(IREFSM,1))

      CALL mma_allocate(CI_VEC,LBLOCK,Label='CI_VEC')
      CALL mma_allocate(SIGMA_VEC,LBLOCK,Label='SIGMA_VEC')

      END SUBROUTINE LUCIA
