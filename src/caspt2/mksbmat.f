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
      SUBROUTINE MKSBMAT()
      use definitions, only: iwp, wp, Byte, u6
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG, VERBOSE
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_global, only: DREF, PREF
      use caspt2_global, only: LUSOLV
      use caspt2_module, only: NASHT
      use pt2_guga, only: NG1, NG2, NG3
      IMPLICIT NONE
C Set up S and B matrices for cases 1..13.

      INTEGER(kind=Byte), ALLOCATABLE :: idxG3(:,:)
      real(kind=wp), ALLOCATABLE:: F1(:), F2(:), F3(:), FD(:), FP(:),
     &                                           G3(:)
      INTEGER(kind=iwp) iLUID
      Logical(kind=iwp) Single_set_of_PCO

      IF(IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,*)' Construct S and B matrices'
      END IF


!     Single_set_of_PCO=.TRUE.
      Single_set_of_PCO=.FALSE.
      If (Single_set_of_PCO) THEN
         CALL MKSMAT()
         CALL MKBMAT()
         RETURN
      END IF

      IF (NASHT/=0) THEN

         CALL mma_allocate(F1,NG1,Label='F1')
         CALL PT2_GET(NG1,'DELTA1',F1)

         CALL mma_allocate(FD,SIZE(DREF),Label='FD')
         CALL MKDREF_RPT2(NASHT,F1,FD,SIZE(DREF))

         CALL mma_deallocate(F1)

         CALL mma_allocate(F2,NG2,Label='F2')
         CALL PT2_GET(NG2,'DELTA2',F2)

         CALL mma_allocate(FP,SIZE(PREF),Label='FP')
         CALL MKPREF_RPT2(NASHT,F2,FP,SIZE(PREF))

         CALL mma_deallocate(F2)

         CALL mma_allocate(F3,NG3,Label='F3')
         CALL PT2_GET(NG3,'DELTA3',F3)

C For the cases A and C, begin by reading in the local storage
C  part of the three-electron density matrix G3:
         CALL mma_allocate(G3,NG3,Label='G3')
         CALL PT2_GET(NG3,'GAMMA3',G3)


         IF (IPRGLB.GE.DEBUG) THEN
           WRITE(u6,'("DEBUG> ",A)') 'CASE SYM S/B-MATRIX NORM'
           WRITE(u6,'("DEBUG> ",A)') '==== === ==============='
         END IF

         CALL mma_allocate(idxG3,6,NG3,label='idxG3')
         iLUID=0
         CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

*
         CALL MKSA(DREF,SIZE(DREF),PREF,SIZE(PREF),NG3,G3,idxG3)
         CALL MKBA(DREF,SIZE(DREF),PREF,SIZE(PREF),FD,FP,NG3,F3,idxG3)

         CALL MKSC(DREF,SIZE(DREF),PREF,SIZE(PREF),NG3,G3,idxG3)
         CALL MKBC(DREF,SIZE(DREF),PREF,SIZE(PREF),FD,FP,NG3,F3,idxG3)

         CALL mma_deallocate(G3)
         CALL mma_deallocate(F3)
         CALL mma_deallocate(idxG3)

         CALL MKSB(DREF,SIZE(DREF),PREF,SIZE(PREF))
         CALL MKBB(DREF,SIZE(DREF),PREF,SIZE(PREF),FD,FP)

         CALL MKSD(DREF,SIZE(DREF),PREF,SIZE(PREF))
         CALL MKBD(DREF,SIZE(DREF),PREF,SIZE(PREF),FD,FP)

         CALL MKSE(DREF,SIZE(DREF))
         CALL MKBE(DREF,SIZE(DREF),FD)

         CALL MKSF(PREF,SIZE(PREF))
         CALL MKBF(DREF,SIZE(DREF),PREF,SIZE(PREF),FP)

         CALL MKSG(DREF,SIZE(DREF))
         CALL MKBG(DREF,SIZE(DREF),FD)

         CALL mma_deallocate(FP)
         CALL mma_deallocate(FD)

      END IF

      CALL MKSH()
      CALL MKBH()

      END SUBROUTINE MKSBMAT
