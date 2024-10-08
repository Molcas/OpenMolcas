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
      SUBROUTINE PREPSD(WFTP,SGS,CIS,LSYM,
     &                  ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,
     &                  NCONF,CI,DET,detocc,detcoeff,SPTRA)
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: SGStruct, CIStruct
      IMPLICIT NONE
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      INTEGER ICNFTAB(*),ISPNTAB(*),ISSTAB(*),IFSBTAB(*)
      INTEGER LSYM,NCONF
      REAL*8 CI(*),DET(*)
      INTEGER IMODE
      CHARACTER(LEN=8) WFTP
      character(len=*), intent(out) :: detocc(NCONF)
      real(8), intent(out) :: detcoeff(NCONF)
      real(8), intent(in) :: SPTRA(*)
      Real*8, Allocatable:: CTMP(:)

C Purpose: Given a RASSCF wave function in Split-GUGA format
C and an orbital transformation matrix for the purpose of
C getting biorthonormal orbitals, prepare a wave function
C in the general SD format, using transformed orbitals.

      IF(WFTP.EQ.'GENERAL ') THEN
C Transform SGUGA to SymmG:
        CALL mma_allocate(CTMP,NCONF,Label='CTMP')
        IMODE=1
        CALL SYG2SGU(IMODE,SGS,CIS,LSYM,ICNFTAB,ISPNTAB,CI,CTMP)
C Transform SymmG to Slater Dets:
        CALL SYGTOSD(ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,CTMP,DET,
     &               detocc,detcoeff,SPTRA)
        CALL mma_deallocate(CTMP)
      ELSE
        DET(1)=CI(1)
      END IF

      END SUBROUTINE PREPSD
