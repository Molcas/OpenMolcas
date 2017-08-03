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

      Real*8 Function x_Golden(ax,bx,cx,f,tol_x,tol_f,xmin,
     &                 rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,
     &                 lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,
     &                 iPrint_Errors)
      Implicit None
      Real*8 :: ax, bx, cx, tol_x, tol_f, xmin
      Real*8, Parameter :: Ratio = 0.5D0*(3.0D0-Sqrt(5.0D0))
      Real*8, Parameter :: RM = 1.0D0-Ratio
      Real*8 :: x1, x2, x3, x4, f2, f3
c External function f and its arguments
      Real*8, External :: f
      Integer :: nij,lMax,nElem
      Real*8 :: rMP(nij,0:nElem),xrMP(nij,nElem),xxrMP(nij,nElem),
     &          xnrMP(nij,nElem),EC(3,nij),AC(3,nij),R_ij(3),C_o_C(3),
     &          Scratch_New(nij*(2+lMax+1)),Scratch_Org(nij*(2+lMax+1))
      Integer :: ij,l,nAtoms,nPert,iPrint_Errors

      x1 = ax
      x4 = cx
      If (Abs(cx-bx) .gt. Abs(bx-ax)) Then
        x2 = bx
        x3 = RM*bx + Ratio*cx
      Else
        x2 = RM*bx + Ratio*ax
        x3 = bx
      End If
      f2 = f(x2,
     &       rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,
     &       nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
      f3 = f(x3,
     &       rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,
     &       nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
      Do While ((Abs(x4-x1) .gt. tol_x*(Abs(x1)+Abs(x2))) .and.
     &          (Abs(f3-f2) .gt. tol_f*(Abs(f2)+Abs(f3))))
        If (f2 .lt. f3) Then
          x4 = x3
          x3 = x2
          f3 = f2
          x2 = RM*x3 + Ratio*x1
          f2 = f(x2,
     &       rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,
     &       nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
        Else
          x1 = x2
          x2 = x3
          f2 = f3
          x3 = RM*x2 + Ratio*x4
          f3 = f(x3,
     &       rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,
     &       nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
        End If
      End Do
      If (f2 .lt. f3) Then
        xmin = x2
        x_Golden = f2
      Else
        xmin = x3
        x_Golden = f3
      End If

      End Function x_Golden
