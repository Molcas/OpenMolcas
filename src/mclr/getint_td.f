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
      SUBROUTINE GETINT_td(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IKSM,JLSM,ICTL,ieaw )
      use Arrays, only: pInt2, KINT2, KINT2a
      use MCLR_Data, only: Square
      use orbinp_mclr, only: NOBPTS
*
* Outer routine for accessing integral block
*
      IMPLICIT None
      Real*8 XINT(*)
      Integer ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IKSM,JLSM,ICTL,ieaw
*
#include "Input.fh"
#include "csm.fh"
       Integer nTest,iXChng,iCoul,nI,nK,nIK,nJ,nL,nJL,nIJ,nKL
*
       NTEST=0
*
*          Write(*,*)'square in getint_td',square
          IF(.not.square) THEN

           IXCHNG=0
           ICOUL=0
           If (ictl.eq.2) ICOUL=1
           If (ictl.eq.3) ICOUL=1
           If (ictl.eq.4) IXCHNG=1
           If (ieaw.ne.0) Then
              CALL GETINC_ABT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                        IXCHNG,IKSM,JLSM,KINT2a,
     &                        pINT2,NSMOB,ICOUL,ieaw )
           Else
              CALL GETINC_ABT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                        IXCHNG,IKSM,JLSM,KINT2,
     &                        pINT2,NSMOB,ICOUL,ieaw )
           End If
          ELSE
           CALL GETINC_ABS_td(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IKSM,JLSM,KINT2,
     &                  pINT2,NSMOB,ICTL)
C
          End If
*
      IF(NTEST.NE.0) THEN
        NI = NOBPTS(ITP,ISM)
        NK = NOBPTS(KTP,KSM)
        IF(IKSM.EQ.0) THEN
          NIK = NI * NK
        ELSE
          NIK = NI*(NI+1)/2
        END IF
        NJ = NOBPTS(JTP,JSM)
        NL = NOBPTS(LTP,LSM)
        IF(JLSM.EQ.0) THEN
          NJL = NJ * NL
        ELSE
          NJL = NJ*(NJ+1)/2
        END IF
        IF(ICOUL.EQ.0) THEN
          WRITE(6,*) ' 2 electron integral block for TS blocks '
          WRITE(6,*) ' Ixchng :', IXCHNG
          WRITE(6,'(1X,4(A,I2,A,I2,A))')
     &    '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &    '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
           CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
        ELSE
          WRITE(6,*) ' Integrals in Coulomb form '
          WRITE(6,'(1X,4(A,I2,A,I2,A))')
     &   '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &   '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
          NIJ = NI*NJ
          NKL = NK*NL
          CALL WRTMAT(XINT,NIJ,NKL,NIJ,NKL)
        END IF

      END IF
*
C     STOP ' Jeppe forced me to stop in GETINT '
      END SUBROUTINE GETINT_td
