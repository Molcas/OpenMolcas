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
      SUBROUTINE PRGASTAB(REST)
      IMPLICIT NONE
      INTEGER REST(*)

      INTEGER IGAS,NGAS
      INTEGER NSYM,ISYM,KORB,KREST

C Executable statements
      WRITE(6,*)
      WRITE(6,*)' GAS restriction table printout'
      WRITE(6,'(A,I5)')'Table size        NSIZE=',REST(1)
      WRITE(6,'(A,I5)')'Table type ID     ITYPE=',REST(2)
      WRITE(6,'(A,I5)')'Nr of partitions  NGAS=',REST(3)
      WRITE(6,'(A,I5)')'Nr of symm labels NSYM =',REST(4)
      WRITE(6,*)' Orbital partitions:'
      NGAS=REST(3)
      NSYM =REST(4)
      KORB=5
      WRITE(6,'(8X,I5,5X,8I5)')
     &      (REST(KORB+ISYM),ISYM=0,NSYM-1)
      DO IGAS=1,NGAS
       WRITE(6,'(I3,5X,I5,5X,8I5)')IGAS,
     &      (REST(KORB+ISYM+(NSYM+1)*IGAS),ISYM=0,NSYM-1)
      END DO
      WRITE(6,*)' Electron population restrictions:'
      KREST=KORB+(NGAS+1)*(NSYM+1)
      WRITE(6,'(5X,A7,5X,30I3)')'Minimum',
     &      (REST(KREST+0+2*(IGAS-1)),IGAS=1,NGAS)
      WRITE(6,'(5X,A7,5X,30I3)')'Maximum',
     &      (REST(KREST+1+2*(IGAS-1)),IGAS=1,NGAS)

      END SUBROUTINE PRGASTAB
