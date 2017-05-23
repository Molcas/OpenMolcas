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
      SUBROUTINE GETINT_MCLR(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,ICOUL ,ieaw)

*
* Outer routine for accessing integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "detdim.fh"
*./CRUN/ : INTIMP used ....
* and NOHSOO (no spin-other-orbit) added by Merethe 19/10-95
#include "crun_mclr.fh"
*./ORBINP/  : NOBPTS used
*.Memory
#include "WrkSpc.fh"
*

#include "Input.fh"
#include "orbinp_mclr.fh"
#include "csm.fh"
#include "genop.fh"
*. Type of operator in action


#include "glbbas_mclr.fh"
      Dimension XINT(*)
*
C      CALL QENTER('GETIN  ')
       NTEST=0
*
          IF(.not.square) THEN
           ip=kint2
           If (ieaw.ne.0) ip=kint2a
           CALL GETINC_ABT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,wORK(ip),
     &                  iWORK(KPINT2),NSMOB,ICOUL,ieaw )
          ELSE
           CALL GETINC_ABS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,wORK(KINT2),
     &                  iWORK(KPINT2),NSMOB,ICOUL )
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
          WRITE(6,'(1H ,4(A,I2,A,I2,A))')
     &    '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &    '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
           CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
        ELSE
          WRITE(6,*) ' Integrals in Coulomb form '
          WRITE(6,'(1H ,4(A,I2,A,I2,A))')
     &   '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &   '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
          NIJ = NI*NJ
          NKL = NK*NL
          CALL WRTMAT(XINT,NIJ,NKL,NIJ,NKL)
        END IF

      END IF
*
C     CALL QEXIT('GETIN ')
      RETURN
      END
