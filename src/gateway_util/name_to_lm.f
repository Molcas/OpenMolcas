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
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine Name_to_lm(BName,l,m)
      Implicit None
      Character(Len=*), Intent(In) :: BName
      Integer, Intent(Out) :: l, m
      Character :: Letter
      Integer :: i, lx, ly, lz
#include "itmax.fh"
#include "angtp.fh"
*
      Letter = BName(2:2)
      Call LoCase(Letter)
      l = 0
      m = 0
      If (Letter .eq. 's') Then
*     Default is s
        Return
      Else if (Letter .eq. 'p') Then
*       p always appear as px, py, pz
        l = 1
        Letter = BName(3:3)
        Call LoCase(Letter)
        If (Letter .eq. 'x') Then
          m = 1
        Else If (Letter .eq. 'y') Then
          m = -1
        Else If (Letter .eq. 'z') Then
          m = 0
        End If
      Else
        l = -1
*       Find if there is an angular label
        Do i=Sum(LBound(AngTp)),Sum(UBound(AngTp))
          If (Letter .eq. AngTp(i)) Then
            l = i
            Exit
          End If
        End Do
        If (l .ge. 0) Then
*         If a label is found it is a spherical shell, just read m
          Read(BName(3:4),*) m
          If (BName(5:5) .eq. '-') m = -m
        Else
*         If no label, this is a Cartesian shell, return -l and some convention for m.
*         We use m=T(ly+lz)-(lx+ly), where T(n) is the nth triangular number: n*(n+1)/2).
*         From here, ly+lz can be recovered as the (int) triangular root of m+l: (sqrt(8*(m+l)+1)-1)/2,
*         lz is m+l-T(ly+lz) and lx is l-(ly+lz). This has the property that all
*         possible combinations of lx,ly,lz are encoded in -l plus a number from -l to l*(l+1)/2,
*         in descending order with priority lx>ly>lz: (3,0,0), (2,1,0), (2,0,1), (1,2,0),
*         (1,1,1), (1,0,2), (0,3,0), (0,2,1), (0,1,2), (0,0,3)
          Read(BName(2:3),*) lx
          Read(BName(4:5),*) ly
          Read(BName(6:7),*) lz
          l = -lx-ly-lz
          m = (ly+lz)*(ly+lz+1)/2-(lx+ly)
        End If
      End If
      Return
*
      End Subroutine Name_to_lm
