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
      FUNCTION NSXFSM(NSMOB,MXPOBS,NO1PS,NO2PS,ISXSM,ADSXA,
     &ISYM,IPRNT)
*
* Number of single excitations of symmetry ISXSM
*
* ISYM = 0 : All symmetry allowed excitations
* ISYM = 1 : Only excitations a+iaj with I.ge.J
* ISYM =-1 : Only excitations a+iaj with I.gt.J
      INTEGER ADSXA(MXPOBS,2*MXPOBS)
      INTEGER NO1PS(*),NO2PS(*)
*
      MSXFSM = 0
C?    WRITE(6,*) ' NSMOB ',NSMOB
      DO 100 IO1SM = 1,NSMOB
        IO2SM = ADSXA(IO1SM,ISXSM)
C?      WRITE(6,*) ' IO1SM,IO2SM',IO1SM,IO2SM
        IF(ISYM.EQ.0.OR.IO1SM.GT.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*NO2PS(IO2SM)
        ELSE IF( ISYM.EQ. 1 .AND. IO1SM.EQ.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*(NO1PS(IO1SM)+1)/2
        ELSE IF( ISYM.EQ.-1 .AND. IO1SM.EQ.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*(NO1PS(IO1SM)-1)/2
        END IF
  100 CONTINUE
*
      NSXFSM = MSXFSM
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)

      IF(NTEST.NE.0) THEN
        WRITE(6,*)
     &  ' Number of single excitations of symmetry ',ISXSM,',',NSXFSM
      END IF
*
      RETURN
      END
* Output
