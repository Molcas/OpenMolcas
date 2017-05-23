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
      SUBROUTINE VSUB(VEC1,IST1,VEC2,IST2,VEC3,IST3,NS)
      implicit none
      integer IST1,IST2,IST3,NS, I, IS1,IS2,IS3
      REAL*8 VEC1,VEC2,VEC3
      DIMENSION VEC1(*),VEC2(*),VEC3(*)
c
      IF(IST1.EQ.1.AND.IST2.EQ.1.AND.IST3.EQ.1)THEN
         DO I=1,NS
            VEC3(I  )=VEC2(I  )-VEC1(I  )
         enddo
      ELSE
         IS1=1
         IS2=1
         IS3=1
         DO I=1,NS
            VEC3(IS3)=VEC2(IS2)-VEC1(IS1)
            IS1=IS1+IST1
            IS3=IS3+IST3
            IS2=IS2+IST2
      enddo
      ENDIF
      RETURN
      END
