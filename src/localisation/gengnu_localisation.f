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
      SubRoutine GenGnu_Localisation(FilNam,Den,Coord,n)
      Implicit None
      Character*12 FilNam
      Integer n
      Real*8  Den(n,n), Coord(3,n)

      Integer  isFreeUnit
      External isFreeUnit

      Character*4  Postfix
      Character*16 FName
      Parameter (Postfix = '.dat')

      Integer i, j, Lu
      Real*8  Dist, D, Fac

C     Open data file.
C     ---------------

      Lu = isFreeUnit(11)
      Write(FName,'(A12,A4)') FilNam,Postfix
      If (FName(1:1) .eq. ' ') FName(1:1) = 'G'
      Do i = 2,12
         If (FName(i:i) .eq. ' ') FName(i:i) = '_'
      End Do
      Call Molcas_Open(Lu,FName)
      Rewind(Lu)

C     Write data file.
C     ----------------

      Fac = 1.0d0/log(1.0d1)
      Do j = 1,n
         Do i = j,n
            Dist = sqrt((Coord(1,j)-Coord(1,i))**2
     &                 +(Coord(2,j)-Coord(2,i))**2
     &                 +(Coord(3,j)-Coord(3,i))**2
     &                 )
            If (abs(Den(i,j)) .lt. 1.0d-16) Then
               D = -4.0d1
            Else
               D = Fac*log(abs(Den(i,j)))
            End If
            Write(Lu,'(1X,1P,D20.10,1X,D20.10)') Dist,D
         End Do
      End Do

C     Close data file.
C     ----------------

      Close(Lu,Status='Keep')

      End
