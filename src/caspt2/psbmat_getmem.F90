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
!
! WRAPPERS FOR PARALLEL S AND B MATRIX ROUTINES
!
      SUBROUTINE PSBMAT_GETMEM(cNAME,lg_M,nSize)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_ga, only: GA_arrays, Allocate_GA_Array
      use constants, only: Zero
      use definitions, only: iwp
      IMPLICIT None
!SVC2010: create square global array S/B for symmetry iSYM
! with integer handle lg_M or if replicate or serial, create
! tridiagonal local array at Work(lg_M)
      Integer(kind=iwp) lg_M, nSize
      CHARACTER(len=*) cNAME

      Integer(kind=iwp) nTri

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_CREATE_STRIPED ('H',nSize,nSize,cNAME,LG_M)
        CALL GA_ZERO (LG_M)
      ELSE
#endif
        nTri=(nSize*(nSize+1))/2
        lg_M=Allocate_GA_Array(nTri,cName)
        GA_Arrays(lg_M)%A(:)=Zero
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE PSBMAT_GETMEM
