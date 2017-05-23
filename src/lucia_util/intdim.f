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
      SUBROUTINE INTDIM(IPRNT)
*
* Number of integrals and storage mode
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* =====
*.Input
* =====
*
#include "mxpdim.fh"
#include "lucinp.fh"
#include "orbinp.fh"
#include "csm.fh"
#include "crun.fh"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
* =======
*. Output
* =======
*
#include "cintfo.fh"
*.1 : Number of one-electron integrals
      NINT1 =  NSXFSM(NSMOB,MXPOBS,NTOOBS,NTOOBS,ITSSX,ADSXA,1,IPRNT)
*.2 : Number of two-electron integrals
        IF(PNTGRP.EQ.1) THEN
*. Full eightfold symmetry can be used
          I12S = 1
          I34S = 1
          I1234S = 1
c        ELSE
*. Only symmetry between 12 and 34
c          I12S = 0
c          I34S = 0
c          I1234S = 1
        END IF
        NINT2 = NDXFSM(NSMOB,NSMSX,MXPOBS,NTOOBS,NTOOBS,NTOOBS,
     &                  NTOOBS,ITSDX,ADSXA,SXDXSX,I12S,I34S,I1234S,
     &                  IPRNT )
*. Number of integrals without complex conjugation symmetry
*. ( used for T1 transformed Hamiltonian)
       I12 = 0
       I34 = 0
       I1234 = 1
       NINT2_NO_CCSYM = NDXFSM(NSMOB,NSMSX,MXPOBS,NTOOBS,NTOOBS,NTOOBS,
     &                  NTOOBS,ITSDX,ADSXA,SXDXSX,I12,I34,I1234,
     &                  IPRNT )
c       IF(ISIMTRH.EQ.1) THEN
c         WRITE(6,*)
c     &   ' Number of two-electron integrals in exp(-T1)Hexp(T1) ',
c     &     NINT2_NO_CCSYM
c       END IF
*. Number of symmetry blocks of one- and two-electron integrals
      NBINT1 = NSMOB
      NBINT2 = NSMOB ** 3
      RETURN
      END
