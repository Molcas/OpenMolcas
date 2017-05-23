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
      SUBROUTINE HSHGET(KEY,KEYDIM,NCOMP,ITEM,
     &                                     NSIZE,ITAB,ITEMID)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ITAB(NSIZE,2)
      DIMENSION KEY(KEYDIM),ITEM(NCOMP,*)
C These parameters determine the hash function:
      PARAMETER (MULT=37,NHASH=997)

      IF(NSIZE.LT.NHASH) GOTO 901
      NULL=ITAB(NSIZE,1)
C Evaluate hash index:
      IND=MOD(KEY(1),NHASH)
      DO I=2,KEYDIM
        IND=MOD(KEY(I)+MULT*IND,NHASH)
      END DO
      IND=IND+1
C IND is a hashed index in interval 1..NHASH < NSIZE
C Find the item with this key:

      LOOKAT=IND
  10  CONTINUE
C Are there (more) items with that hash signature?
      IF(ITAB(LOOKAT,1).EQ.NULL) GOTO 30
C Try to identify an item which has the given key:
      ITEMID=ITAB(LOOKAT,2)
      DO I=1,KEYDIM
        IF(ITEM(I,ITEMID).NE.KEY(I)) GOTO 20
      END DO
C Here, if we have identified the item.
      RETURN

C Here, if we have not yet identified the item.
  20  CONTINUE
      LOOKAT=ITAB(LOOKAT,1)
      GOTO 10

C Here, if we have failed to find such an item.
  30  CONTINUE
      ITEMID=0
      RETURN
 901  CONTINUE
      WRITE(6,*)' HSHGET: Table size must be at least as'
      WRITE(6,*)'         big as NHASH, presently =', NHASH
      CALL ABEND()
      END
      SUBROUTINE HSHPUT(KEYDIM,NCOMP,ITEM,NSIZE,ITAB,ITEMID)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ITAB(NSIZE,2)
      DIMENSION ITEM(NCOMP,*)
C These parameters determine the hash function:
      PARAMETER (MULT=37,NHASH=997)

      IF(NSIZE.LT.NHASH) GOTO 901
      NULL=ITAB(NSIZE,1)
      IFREE=ITAB(NSIZE,2)
CTEST      write(*,*)' HSHPUT: NSIZE=',NSIZE
CTEST      write(*,*)' HSHPUT: NULL =',NULL
CTEST      write(*,*)' HSHPUT: IFREE=',IFREE
      IF(ITAB(IFREE,1).EQ.NULL) GOTO 902
C Evaluate hash index:
      IND=MOD(ITEM(1,ITEMID),NHASH)
      DO I=2,KEYDIM
        IND=MOD(ITEM(I,ITEMID)+MULT*IND,NHASH)
      END DO
      IND=IND+1
CTEST      write(*,*)' HSHPUT:   IND=',  IND
C IND is a hashed index in interval 1..NHASH < NSIZE
C Find the last item with this key:

      LOOKAT=IND
  10  CONTINUE
C Are there already items with that hash signature?
CTEST      write(*,*)' HSHPUT:LOOKAT=',LOOKAT
      IF(ITAB(LOOKAT,1).EQ.NULL) GOTO 30
CTEST      write(*,*)'   Already taken.'
      LOOKAT=ITAB(LOOKAT,1)
      GOTO 10

  30  CONTINUE
C No more items with the same signature.
C Put the new item in the table at a free location.
CTEST      write(*,'(1x,A,2I8)')
CTEST     &       ' ITAB(LOOKAT..=',ITAB(LOOKAT,1),ITAB(LOOKAT,2)
CTEST      write(*,*)' Set 1st pointer to the IFREE location:'
      ITAB(LOOKAT,1)=IFREE
      ITAB(LOOKAT,2)=ITEMID
      NEXT=ITAB(IFREE,1)
      ITAB(IFREE,1)=NULL
      ITAB(NSIZE,2)=NEXT
CTEST      write(*,*)' The first few free ITAB elements:'
c      DO I=IFREE-3,IFREE+3
CTEST        write(*,'(1x,3I8)') I,ITAB(I,1),ITAB(I,2)
c      END DO
      RETURN
 901  CONTINUE
      WRITE(6,*)' HSHPUT: Table size must be at least as'
      WRITE(6,*)'         big as NHASH, presently =', NHASH
      CALL ABEND()
 902  CONTINUE
      WRITE(6,*)' HSHPUT: Table is already full.'
      WRITE(6,*)' Size NSIZE is too small, NSIZE =', NSIZE
      CALL ABEND()
      END
      SUBROUTINE HSHINI(NSIZE,ITAB,NULL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ITAB(NSIZE,2)
C These parameters determine the hash function
      PARAMETER (NHASH=997)

      IF(NSIZE.LT.NHASH) GOTO 901
      DO I=1,NHASH
        ITAB(I,1)=NULL
        ITAB(I,2)=NULL
      END DO
      IFREE=NHASH+1
      DO I=IFREE,NSIZE-1
        ITAB(I,1)=I+1
        ITAB(I,2)=NULL
      END DO
      ITAB(NSIZE,1)=NULL
      ITAB(NSIZE,2)=IFREE
CTEST      write(*,*)' The first few free ITAB elements:'
CTEST      DO I=IFREE-3,IFREE+3
CTEST        write(*,'(1x,3I8)') I,ITAB(I,1),ITAB(I,2)
CTEST      END DO
CTEST      write(*,*)' The last few ITAB elements:'
CTEST      DO I=NSIZE-5,NSIZE
CTEST        write(*,'(1x,3I8)') I,ITAB(I,1),ITAB(I,2)
CTEST      END DO
      RETURN
 901  CONTINUE
      WRITE(6,*)' HSHINI: Table size must be at least as'
      WRITE(6,*)'         big as NHASH, presently =', NHASH
      CALL ABEND()
      END
