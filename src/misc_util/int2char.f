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
* Copyright (C) 2003, Oleh Danyliv                                     *
************************************************************************
      CHARACTER*(*) FUNCTION INT2CHAR(N,Base)
************************************************************************
*                                                                      *
* Object: to convert integer into string                               *
*                                                                      *
* Called from: rdctl                                                   *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Oleh Danyliv, 2003                                       *
*                                                                      *
************************************************************************
c     Converts integer number N into character length of Base
c      Input:
c       N - integer for conversion
c       Base - length of character (<=100)
c      Output
c       INT2CHAR - converted i
      IMPLICIT NONE
      INTEGER N, Base, N_h
      INTEGER i,j,i_h
      INTEGER Shift
      CHARACTER*100 str
      CHARACTER*1 str_m(100)
      EQUIVALENCE (str,str_m)
      LOGICAL First0

      Shift = 48
      N_h = N
      i_h = 0
      First0 = .TRUE.
      str = ' '
      DO i=1, Base
         j  = INT( N_h/(10.d0**(Base-i)) )
         IF(First0.AND.(j.EQ.0)) THEN
            i_h = i_h + 1
         ELSE
            First0 = .FALSE.
            str_m(i-i_h) = CHAR(j+Shift)
         END IF
         N_h = N_h - j*(10**(Base-i))
      END DO


      INT2CHAR = str

      RETURN
      END
