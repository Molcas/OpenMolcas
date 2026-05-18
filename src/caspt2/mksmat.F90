!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
      SUBROUTINE MKSMAT()
      use definitions, only: iwp, wp, u6, byte
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG, VERBOSE
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_global, only: DREF, PREF, LUSOLV
      use caspt2_module, only: NASHT
      use caspt2_module, only: NG3
      IMPLICIT None
!     Set up S matrices for cases 1..13.
      INTEGER(kind=byte), ALLOCATABLE :: idxG3(:,:)

      real(kind=wp), ALLOCATABLE:: G3(:)
      integer(kind=iwp) iLUID


      IF(IPRGLB.GE.VERBOSE) THEN
        WRITE(u6,*)
        WRITE(u6,*)' Construct S matrices'
      END IF

      IF(NASHT.GT.0) THEN
!SVC: print header for debug info
        IF(IPRGLB.GE.DEBUG) THEN
          WRITE(u6,'("DEBUG> ",A)') 'CASE SYM S-MATRIX NORM'
          WRITE(u6,'("DEBUG> ",A)') '==== === ============='
        END IF
! For the cases A and C, begin by reading in the local storage
!  part of the three-electron density matrix G3:
        CALL mma_allocate(G3,NG3,Label='G3')
        CALL PT2_GET(NG3,'GAMMA3',G3)

        CALL mma_allocate(idxG3,6,NG3,label='idxG3')
        iLUID=0
        CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

        CALL MKSA(DREF,SIZE(DREF),PREF,SIZE(PREF),NG3,G3,idxG3)
        CALL MKSC(DREF,SIZE(DREF),PREF,SIZE(PREF),NG3,G3,idxG3)

        CALL mma_deallocate(G3)
        CALL mma_deallocate(idxG3)

!-SVC20100902: For the remaining cases that do not need G3, use replicate arrays
        CALL MKSB(DREF,SIZE(DREF),PREF,SIZE(PREF))
        CALL MKSD(DREF,SIZE(DREF),PREF,SIZE(PREF))
        CALL MKSE(DREF,SIZE(DREF))
        CALL MKSF(PREF,SIZE(PREF))
        CALL MKSG(DREF,SIZE(DREF))
      END IF

      Call MKSH()

      END SUBROUTINE MKSMAT
