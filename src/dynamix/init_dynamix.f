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
C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
      SUBROUTINE Init_Dynamix
      IMPLICIT REAL*8 (a-h,o-z)
#include "MD.fh"
#include "WrkSpc.fh"
C      CHARACTER*180 Lines(10)
*
C#ifdef _SKIP_
C      Lines(1 )=_MOLCAS_VERSION_
C#ifdef _DEMO_
C      Lines(2 )='DEMO VERSION '
C#else
C      Lines(2 )=' '
C#endif
c      Lines(3 )=' '
C      Lines(4 )='D Y N A M I X'
C      Lines(5 )=' '
C      Lines(6 )='Written by Igor Schapiro'
C      Lines(7 )=' '
C      Lines(8 )='A program for molecular dynamics simulations'
C      Lines(9 )=' '
C      Lines(10)='Code compiled :'//_BUILD_DATE_
C      lLine=Len(Lines(1))
c      CALL Banner(Lines,10,lLine)
*
C      WRITE (6,*)
C      WRITE (6,'(A,I8,A)')   ' The Dynamix program has ',
C     &      mxMem,' double precision words memory available.'
C      WRITE (6,'(A,I8,A)')   '                       ',
C     &      mxMem*8/(1024**2),' MB'
C      WRITE (6,*)
C#endif
C
C     Set the defalut values
C
      THERMO=0
      TEMP=298.15
      VELO=0
      POUT=0
      iPrint=2
      DT=1.0D1
      TIME=0.0D0
      RESTART=0.0D0
      lH5Restart= .False.
*
*      Etot0=0.0D0
*
      RETURN
*
      END
