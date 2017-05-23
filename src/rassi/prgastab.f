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
      SUBROUTINE PRGASTAB(LREST)
      IMPLICIT NONE
      INTEGER IGAS,LREST,NGAS
      INTEGER NSYM,ISYM,KORB,KREST
#include "WrkSpc.fh"

C Executable statements
      WRITE(6,*)
      WRITE(6,*)' GAS restriction table printout'
      WRITE(6,'(A,I5)')'Table size        NSIZE=',IWORK(LREST+0)
      WRITE(6,'(A,I5)')'Table type ID     ITYPE=',IWORK(LREST+1)
      WRITE(6,'(A,I5)')'Nr of partitions  NGAS=',IWORK(LREST+2)
      WRITE(6,'(A,I5)')'Nr of symm labels NSYM =',IWORK(LREST+3)
      WRITE(6,*)' Orbital partitions:'
      NGAS=IWORK(LREST+2)
      NSYM =IWORK(LREST+3)
      KORB=5
      WRITE(6,'(8X,I5,5X,8I5)')
     &      (IWORK(LREST-1+KORB+ISYM),ISYM=0,NSYM)
      DO IGAS=1,NGAS
       WRITE(6,'(I3,5X,I5,5X,8I5)')IGAS,
     &      (IWORK(LREST-1+KORB+ISYM+(NSYM+1)*IGAS),ISYM=0,NSYM)
      END DO
      WRITE(6,*)' Electron population restrictions:'
      KREST=KORB+(NGAS+1)*(NSYM+1)
      WRITE(6,'(5X,A7,5X,30I3)')'Minimum',
     &      (IWORK(LREST-1+KREST+0+2*(IGAS-1)),IGAS=1,NGAS)
      WRITE(6,'(5X,A7,5X,30I3)')'Maximum',
     &      (IWORK(LREST-1+KREST+1+2*(IGAS-1)),IGAS=1,NGAS)
      RETURN
      END
