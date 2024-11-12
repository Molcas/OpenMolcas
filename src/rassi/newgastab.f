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
      SUBROUTINE NEWGASTAB(NSYM,NGAS,NGASORB,NGASLIM,ICASE)
      use stdalloc, only: mma_allocate
      use rassi_global_arrays, only: REST1, REST2, REST
      IMPLICIT NONE
      INTEGER NSYM,NGAS,NGASORB(NSYM,NGAS),NGASLIM(2,NGAS), ICASE

      INTEGER NSIZE,ITYPE
      INTEGER IGAS,ISUM,ISYM,KORB,KREST,LPOS

C Executable statements
      NSIZE=4+(NGAS+1)*(NSYM+1)+2*NGAS
      ITYPE=91
      SELECT CASE (ICASE)
      CASE (1)
         CALL mma_allocate(REST1,NSIZE,Label='REST1')
         REST=>REST1(:)
      CASE (2)
         CALL mma_allocate(REST2,NSIZE,Label='REST2')
         REST=>REST2(:)
      CASE DEFAULT
         Write(6,*) 'NEWGASTAB: Illegal ICASE value'
         Write(6,*) 'ICASE=',ICASE
      END SELECT
      REST(1)=NSIZE
      REST(2)=ITYPE
      REST(3)=NGAS
      REST(4)=NSYM
      KORB=5
CTEST      write(*,*)' In NEWGASTAB. NGASORB array is:'
CTEST      do igas=1,ngas
CTEST        write(*,'(1x,8i5)')(ngasorb(isym,igas),isym=1,nsym)
CTEST      end do
CTEST      write(*,*)' In NEWGASTAB. NGASLIM array is:'
CTEST      write(*,'(1x,20i3)')(ngaslim(1,igas),igas=1,ngas)
CTEST      write(*,'(1x,20i3)')(ngaslim(2,igas),igas=1,ngas)

      DO IGAS=1,NGAS
       ISUM=0
       DO ISYM=1,NSYM
        LPOS=KORB+ISYM+(NSYM+1)*IGAS
        REST(LPOS)=2*NGASORB(ISYM,IGAS)
        ISUM=ISUM+2*NGASORB(ISYM,IGAS)
       END DO
       LPOS=KORB+0+(NSYM+1)*IGAS
       REST(LPOS)=ISUM
      END DO
      DO ISYM=0,NSYM
       ISUM=0
       DO IGAS=1,NGAS
        LPOS=KORB+ISYM+(NSYM+1)*IGAS
        ISUM=ISUM+REST(LPOS)
       END DO
       LPOS=KORB+ISYM
       REST(LPOS)=ISUM
      END DO
      KREST=KORB+(NGAS+1)*(NSYM+1)
      LPOS=KREST
      DO IGAS=1,NGAS
       REST(LPOS  )=NGASLIM(1,IGAS)
       REST(LPOS+1)=NGASLIM(2,IGAS)
       LPOS=LPOS+2
      END DO

      nullify(REST)
      END SUBROUTINE NEWGASTAB
