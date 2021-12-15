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
      Subroutine PrList(TEXT,Char,nDim,FI,N1,N2)
************************************************************************
*                                                                      *
*     Object: To generate a cartesian output with atomic labels        *
*             N1 and N2 are the real limits of dummy FI                *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Character*(*) TEXT
      Character*(*) Char(nDim)
      Real*8 FI(N1,N2)
*
*
      Lu=6
      WRITE (Lu,100) Text
100   FORMAT (//,1X,A,/)
      WRITE (Lu,200)
200   FORMAT (5X,'ATOM',21X,'X',19X,'Y',19X,'Z',/)
      Do 10 I = 1, NDIM
         If (N1.EQ.3) Then
            Write (Lu,300) Char(I),(FI(J,I),J=1,3)
300         Format (5X,A,3X,3F20.10)
         Else
            Write (Lu,300) Char(I),(FI(I,J),J=1,3)
         End If
10    Continue
*
      Return
      End
