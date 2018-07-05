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
* Copyright (C) 2018, Ignacio Fdez. Galvan                             *
************************************************************************
*  Clean_BName
*
*> @brief
*>   Clean a basis function name for printing.
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Clean up a basis function name for printing purposes, by removing
*> unneeded zeros and padding.
*>
*> @param[in] BName  Basis function name
*> @param[in] Offset Initial characters not to be touched
*>
*> @return Cleaned basis function name
************************************************************************
      Function Clean_BName(BName,Offset)
      Implicit None
#include "Molcas.fh"
      Character(Len=*), Intent(In) :: BName
      Integer, Intent(In) :: Offset
      Character(Len=LENIN8) :: Clean_BName
      Character(Len=8) :: Clean
      Integer :: i,Err
      Logical :: Cart
*
      Clean=BName(Offset+1:)
* For spherical functions, the 3rd character is a letter
*   (counting s and p as spherical),
* for Cartesian functions, it is a number
      Read(Clean(3:3),'(I1)',IOStat=Err) i
      Cart=(Err.eq.0)
*
      If (Cart) Then
* For Cartesian functions, remove zeros if all indices are below 10,
* and shift everything one space to the right
        If ((Clean(2:2).eq.'0').And.
     &      (Clean(4:4).eq.'0').And.
     &      (Clean(6:6).eq.'0'))
     &     Clean(2:)=Clean(3:3)//Clean(5:5)//Clean(7:7)
        Clean(2:)=Clean(1:)
        Clean(1:1)=' '
      Else
* For spherical functions, remove leading zero in m number,
* replace leading zero in shell number with space,
* replace "*0" tag with a single *
        If (Clean(1:1).eq.'0') Clean(1:1)=' '
        If (Clean(1:2).eq.'*0') Clean(1:2)=' *'
        If (Clean(4:4).eq.'0') Clean(4:)=Clean(5:)
      End If
*
      Clean_BName=BName(1:Offset)//Clean
      Return
*
      End Function Clean_BName
