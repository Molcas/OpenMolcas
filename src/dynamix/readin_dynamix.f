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
      SUBROUTINE ReadIn_Dynamix(Task,nTasks,mTasks)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Task(nTasks)
*
C
C     Copy input from standard input to a local scratch file
C
      LuSpool=isfreeunit(21)
      CALL SpoolInp(LuSpool)
C
C     Read input
C
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix calls RdInp_Dynamix.'
#endif
      CALL RdInp_Dynamix(LuSpool,Task,nTasks,mTasks)
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix back from RdInp_Dynamix.'
#endif
C
C     Remove local copy of standard input
C
      CLOSE(LuSpool)
*
      RETURN
*
      END
