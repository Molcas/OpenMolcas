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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************

C*MODULE MTHLIB  *DECK PRTRIL
      SUBROUTINE PRTRIL(D,N,mode)
C
      IMPLICIT real*8 (A-H,O-Z)
C
      DIMENSION D(*)
      integer, optional :: mode
C
      if (mode.eq.1) then
       do i = 1, n
         d(i*(i+1)/2) = d(i*(i+1)/2)*2.0d+00
        end do
      end if
      MAX = 5
      MM1 = MAX - 1
      DO 120 I0=1,N,MAX
         IL = MIN(N,I0+MM1)
         WRITE(6,9008)
         WRITE(6,9028) (I,I=I0,IL)
         WRITE(6,9008)
         IL = -1
         DO 100 I=I0,N
            IL=IL+1
            J0=I0+(I*I-I)/2
            JL=J0+MIN(IL,MM1)
            if (mode.eq.1) then
            WRITE(6,9048) I,'        ',(D(J)*0.5d+00,J=J0,JL)
            else
            WRITE(6,9048) I,'        ',(D(J),J=J0,JL)
            end if
  100    CONTINUE
  120 CONTINUE
      if (mode.eq.1) then
       do i = 1, n
         d(i*(i+1)/2) = d(i*(i+1)/2)*0.5d+00
        end do
      end if
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(15X,10(4X,I4,3X))
 9048 FORMAT(I5,2X,A8,10F11.6)
      END
