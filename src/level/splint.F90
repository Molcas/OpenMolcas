!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
      SUBROUTINE SPLINT(LNPT,NTP,R1,V1,MBEG,MEND,XX,YY)
!** Subroutine to generate (if LNPT.ge.0) 4*NTP coefficients CSP(J)
!  of a cubic spline passing through the NTP points (R1(J),V1(J))
!  and to then calculate values of the resulting function YY(I) at the
!  entering abscissae values XX(I) for  I=MBEG to MEND.
!** If LNPT < 0 , generate function values at the given XX(I) using
!  the coefficients CSP(J) obtained and SAVEd on a preceding call.
!** Assumes both R1(J) & XX(I) are monotonic increasing.
!+++++ Calls only subroutines SPLINE and PLYINTRP ++++++++++++++++++++++
!=======================================================================
      INTEGER MAXSP
      PARAMETER (MAXSP=6400)
      INTEGER  I,IER,I1ST,IDER,JK,K,KK,LNPT,N2,N3,NIPT,NTP,MBEG,MEND
      REAL*8 EPS,R2,RI,RRR,TTMP,R1(NTP),V1(NTP),CSP(MAXSP),             &
     &  YY(MEND),XX(MEND)
      SAVE CSP
!
      JK = 0
      IF(4*NTP.GT.MAXSP) THEN
          WRITE(6,602) MAXSP,NTP
!         STOP
          CALL ABEND()
          ENDIF
      EPS= 1.D-6*(R1(2)-R1(1))
      N2= 2*NTP
      N3= 3*NTP
      IF(LNPT.GT.0) THEN
!** On first pass for a given data set, generate spline function
!  coefficients in subroutine SPLINE
!** Start by using a cubic polynomial at each end of the range to get
!  the first derivative at each end for use in defining the spline.
          IDER= 1
          NIPT= 4
          I1ST= NTP-3
          CALL PLYINTRP(R1(I1ST),V1(I1ST),NIPT,R1(NTP),CSP,NIPT,IDER)
          TTMP= CSP(2)
          CALL PLYINTRP(R1,V1,NIPT,R1(1),CSP,NIPT,IDER)
          CSP(1)= CSP(2)
          CSP(2)= TTMP
!** Now call routine to actually generate spline coefficients
          CALL LEVEL_SPLINE(R1,V1,NTP,3,CSP,MAXSP,IER)
          IF(IER .NE. 0) THEN
              WRITE(6,604)
!             STOP
              CALL ABEND()
              ENDIF
          ENDIF
      IF(MEND.LT.MBEG) GO TO 99
!** Now, use spline to generate function at desired points XX(I)
      DO  I= MBEG,MEND
          RI= XX(I)
          RRR= RI-EPS
          KK= 1
!** For a monotonic increasing distance array XX(I),  this statement
!  speeds up the search for which set of cubic coefficients to use.
          IF(I.GT.MBEG) THEN
              IF(XX(I).GT.XX(I-1)) KK= JK
              ENDIF
          DO  K= KK,NTP
              JK= K
              IF(R1(K).GE.RRR) GO TO 64
              ENDDO
   64     CONTINUE
          JK= JK-1
          IF(JK.LT.1) JK= 1
          R2= RI-R1(JK)
          YY(I)= CSP(JK)+R2*(CSP(NTP+JK)+R2*(CSP(N2+JK)+R2*CSP(N3+JK)))
          ENDDO
   99 RETURN
  602 FORMAT(' *** ERROR in SPLINT ***  Array dimension  MAXSP=',I4,    &
     &  ' cannot contain spline coefficients for  NTP=',I4)
  604 FORMAT(' *** ERROR in generating spline coefficients in SPLINE')
      END
