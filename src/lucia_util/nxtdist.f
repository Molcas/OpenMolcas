************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2012, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE NXTDIST(NSMST,NGRP,NELMNT,KGRP,ISMDFGP,SCR,
     &                   NACTSYM,NONEW)
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
*Giovanni Li Manni Feb 2012
*
**************
*. Input
**************
* NSMST   = Number of Irreps
* NGRP    = Number of Groups
* NELMNT  = Number of GAS spaces
* KGRP    = K supergroup
* ISMDFGP = Symmetry distributions per group
**************
* OUTPUT
**************
* NONEW
********************
*. Input and output
********************
      INTEGER KGRP(NELMNT),ISMDFGP(NSMST,NGRP),NELMNT
      INTEGER SCR(NELMNT), NACTSYM(NGRP)
*
       NTEST = 0
       IF(NTEST.NE.0) THEN
         write(6,*) 'NACTSYM : '
         call IWRTMA(NACTSYM,1,NGRP,1,NGRP)
         write(6,*) ' ISMDFGP  Table:'
         do J = 1, NSMST
           WRITE(6,'(40I2)') (ISMDFGP(J,I), I=1,NGRP)
         end do
         write(6,*) 'Initial SCR'
         call IWRTMA(SCR,1,NELMNT,1,NELMNT)
       END IF
*
      IF(NELMNT.EQ.0) THEN
        NONEW = 1
        GOTO 1001
      END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF(SCR(IPLACE).LT.NACTSYM(KGRP(IPLACE))) THEN
          SCR(IPLACE) =  SCR(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
          DO JPLACE = 1, IPLACE
            SCR(JPLACE) = 1
          END DO
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF(NTEST.ne.0) then
        write(6,*) 'New ISMDFGP'
        write(6,'(40I2)') (ISMDFGP(SCR(IGAS),KGRP(IGAS)), IGAS=1,NELMNT)
      end if

      IF(NTEST.NE.0) THEN
        write(6,*) ' New SCR'
        CALL IWRTMA(SCR,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
