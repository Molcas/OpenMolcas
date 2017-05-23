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
      INTEGER FUNCTION NEWGASTAB(NSYM,NGAS,NGASORB,NGASLIM)
      IMPLICIT NONE
      INTEGER LREST,NSIZE,ITYPE
      INTEGER NSYM,NGAS,NGASORB(NSYM,NGAS),NGASLIM(2,NGAS)
      INTEGER IGAS,ISUM,ISYM,KORB,KREST,LPOS
#include "WrkSpc.fh"

C Executable statements
      NSIZE=4+(NGAS+1)*(NSYM+1)+2*NGAS
      ITYPE=91
      CALL GETMEM('GasTab','Allo','Inte',LREST,NSIZE)
      IWORK(LREST+0)=NSIZE
      IWORK(LREST+1)=ITYPE
      IWORK(LREST+2)=NGAS
      IWORK(LREST+3)=NSYM
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
        LPOS=LREST-1+KORB+ISYM+(NSYM+1)*IGAS
        IWORK(LPOS)=2*NGASORB(ISYM,IGAS)
        ISUM=ISUM+2*NGASORB(ISYM,IGAS)
       END DO
       LPOS=LREST-1+KORB+0+(NSYM+1)*IGAS
       IWORK(LPOS)=ISUM
      END DO
      DO ISYM=0,NSYM
       ISUM=0
       DO IGAS=1,NGAS
        LPOS=LREST-1+KORB+ISYM+(NSYM+1)*IGAS
        ISUM=ISUM+IWORK(LPOS)
       END DO
       LPOS=LREST-1+KORB+ISYM
       IWORK(LPOS)=ISUM
      END DO
      KREST=KORB+(NGAS+1)*(NSYM+1)
      LPOS=LREST-1+KREST
      DO IGAS=1,NGAS
       IWORK(LPOS  )=NGASLIM(1,IGAS)
       IWORK(LPOS+1)=NGASLIM(2,IGAS)
       LPOS=LPOS+2
      END DO

      NEWGASTAB=LREST
      RETURN
      END
