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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************
      SUBROUTINE ORDSTR_MCLR(IINST,IOUTST,NELMNT,ISIGN,IPRNT)
!
! ORDER A STRING OF INTEGERS TO ASCENDING ORDER
!
! IINST : INPUT STRING IS IINST
! IOUTST : OUTPUT STRING IS IOUTST
! NELMNT : NUMBER OF INTEGERS IN STRING
! ISIGN :  SIGN OF PERMUTATION : + 1 : EVEN PERMUTATIONN
!                                - 1 : ODD  PERMUTATION
!
! THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
! ( HE IS HEREBY AKNOWLEDGED , AND I AM EXCUSED )
!
! IMPLEMENTED MORE TRANSPARENT BUBBLE SORTING INSTEAD
!               JR NOV 2006
!
      IMPLICIT None
      Integer NELMNT
      Integer IINST(NELMNT),IOUTST(NELMNT)
      Integer ISIGN,IPRNT

      Integer iTemp, iPass, I

      IF(NELMNT.EQ.0) RETURN

      ISIGN=1
      iTEMP=0

10       iPass=0
       DO I=1,NELMNT-1
       IF(IINST(I).GT.IINST(I+1)) THEN
         iTEMP=IINST(I)
         IINST(I)=IINST(I+1)
         IINST(I+1)=iTEMP
         ISIGN=-1*ISIGN
         iPass=1
         ENDIF
       ENDDO
       IF(IPASS.NE.0) GOTO 10


              DO I=1,NELMNT
              IOUTST(I)=IINST(I)
        ENDDO
!
!      NTEST = 0
!      NTEST = MAX(NTEST,IPRNT)
!
!      CALL iCOPY(NELMNT,IINST,1,IOUTST,1)
!      ISIGN = 1
!
!
!        JOE = 1
!10      I = JOE
!20      CONTINUE
!        IF(I.EQ.NELMNT) GO TO 50
!        IF(IOUTST(I).LE.IOUTST(I+1)) GO TO 40
!        JOE = I + 1
!30      iSWAP = IOUTST(I)
!        ISIGN = - ISIGN
!        IOUTST(I) = IOUTST(I+1)
!        IOUTST(I+1) = iSWAP
!       IF(I.EQ.1) GO TO 10
!        I = I - 1
!        IF(IOUTST(I).GT.IOUTST(I+1)) GO TO 30
!        GO TO 10
!40      I = I + 1
!      GO TO 20
!
!     END ORDER
!
!50    CONTINUE
!      IF( NTEST .GE.200) THEN
!        WRITE(6,*)  ' INPUT STRING ORDERED STRING ISIGN ',NELMNT
!        CALL IWRTMA(IINST,1,NELMNT,1,NELMNT)
!        CALL IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
!        WRITE(6,*) ' ISIGN : ', ISIGN
!      END IF
!
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRNT)
      END SUBROUTINE ORDSTR_MCLR
