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
* Copyright (C) 1995, Bernd Artur Hess                                 *
************************************************************************
C $Id: relsewc.r,v 1.4 1995/05/08 14:08:53 hess Exp $
C calculate relativistic operators
C   Bernd Artur Hess, hess@uni-bonn.de
C
C
      SUBROUTINE SCFCLI4(idbg,epsilon,S,H,REVTA,
     *SINVA,NA,NB,ISIZEA,ISIZEB,VELIT,
     *AAA,AAB,ISYMA,ISYMB,CMM1,CMM2,
     *EV4,BU2,BU6,EIG4,EW4,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(ISIZEA),P(ISIZEA),
     *          EV4(ISIZEA),S(ISIZEA)
      DIMENSION SINVA(NA,NA),
     *          REVTA(NA,NA),
     *          BU2(NA,NB),
     *          AAA(NA),AAB(NB)
      DIMENSION CMM1(NB,NA),CMM2(NA,NB)
      DIMENSION EW4(NA),EIG4(NA,NA),BU6(NA,NA)
C
*     Call qEnter('scfcli4')
C
C
      DO I=1,NB
      DO J=1,NA
      CMM1(I,J)=-CMM2(J,I)
      ENDDO
      ENDDO
C
C
C
C     Now we changed sign due to the imaginary character
C
      IJ=0
      DO I=1,NA
      DO J=1,I
      IJ=IJ+1
      EV4(IJ)=0.0D0
      DO L=1,NB
      EV4(IJ)=EV4(IJ)-CMM2(I,L)*CMM1(L,J)
      ENDDO
      ENDDO
      ENDDO
C
C
*     write(*,*)
*     write(*,*) 'Final BSS matrix no Hess part alpha**2 yet'
*     write(*,*) EV4
C
      DO IJ=1,ISIZEA
      EV4(IJ) = 0.5D0*(1.0D0/(VELIT*VELIT))*EV4(IJ)
      ENDDO
C
      CALL AddMar(ISIZEA,EV4,H)
C
C
C
C
culf
      if(idbg.gt.0)CALL PRMAT(IDBG,h,na,nb,'h   oper')
      CALL Sogr(idbg,NA,S,SINVA,P,BU6,EW4)
C
      CALL Diagr(H,NA,EIG4,EW4,SINVA,BU6,EV4)
C
*     if(idbg.gt.0)WRITE (idbg,*) '--- EIGENVALUES OF H MATRIX ---'
*     if(idbg.gt.0)WRITE (idbg,'(4D20.12)') EW4
C
*     write(*,*) 'END OF SCFCLI4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
C
*     Call qExit('scfcli4')
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_real(epsilon)
        CALL Unused_real_array(REVTA)
        CALL Unused_integer(ISIZEB)
        CALL Unused_real_array(AAA)
        CALL Unused_real_array(AAB)
        CALL Unused_integer(ISYMA)
        CALL Unused_integer(ISYMB)
        CALL Unused_real_array(BU2)
      END IF
      END
