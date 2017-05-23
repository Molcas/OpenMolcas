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
      Subroutine Order_Arrays(mode,A1,N1,N2,P2,SCR1)

      Implicit Real*8 (a-h,o-z)
      Character*4 mode
      Integer N1, N2
      Real*8  A1(N1,N2), P2(N2), SCR1(N1)

      If (mode.eq.'decr') Then
         Do j=1,N2-1
            Do k=j+1,N2
               If (P2(j).lt.P2(k)) Then
                  tmp=P2(j)
                  P2(j)=P2(k)
                  P2(k)=tmp
                  call dcopy_(N1,A1(1,j),1,SCR1(1),1)
                  call dcopy_(N1,A1(1,k),1,A1(1,j),1)
                  call dcopy_(N1,SCR1(1),1,A1(1,k),1)
               EndIf
            End Do
         End Do
      ElseIf (mode.eq.'incr') Then
         Do j=1,N2-1
            Do k=j+1,N2
               If (P2(j).gt.P2(k)) Then
                  tmp=P2(j)
                  P2(j)=P2(k)
                  P2(k)=tmp
                  call dcopy_(N1,A1(1,j),1,SCR1(1),1)
                  call dcopy_(N1,A1(1,k),1,A1(1,j),1)
                  call dcopy_(N1,SCR1(1),1,A1(1,k),1)
               EndIf
            End Do
         End Do
      Else
         write(6,*) ' In routine Order_Arrays: wrong mode! '
         Call Abend
      EndIf

      Return
      End
