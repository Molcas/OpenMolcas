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
      Subroutine SqToTri_Q(SqMat,TriMat,iDi)
      Implicit Real*8 (a-h,o-z)
      Dimension SqMat(*),TriMat(*)
      kaunter=0
      Do 10, i=1,iDi
        Do 20, j=1,i
          kaunter=kaunter+1
          TriMat(kaunter)=SqMat(i+(j-1)*iDi)
20      Continue
10    Continue
      Return
      End
