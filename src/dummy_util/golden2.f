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

#ifndef _HAVE_EXTRA_

      Real*8 Function Golden2(ax,bx,cx,f,tol_x,tol_f,xmin,
     &                 q_a,q_b,dipole_a,dipole_b,r_a,r_b)
      Implicit None
      Real*8 :: ax, bx, cx, tol_x, tol_f, xmin
      Real*8, Parameter :: Ratio = 0.5D0*(3.0D0-Sqrt(5.0D0))
      Real*8, Parameter :: RM = 1.0D0-Ratio
      Real*8 :: x1, x2, x3, x4, f2, f3
c External function f and its arguments
      Real*8, External :: f
      Real*8 :: q_a,q_b,dipole_a,dipole_b,r_a,r_b
      Logical, Parameter :: Absolute=.True.

      x1 = ax
      x4 = cx
      If (Abs(cx-bx) .gt. Abs(bx-ax)) Then
        x2 = bx
        x3 = RM*bx + Ratio*cx
      Else
        x2 = RM*bx + Ratio*ax
        x3 = bx
      End If
      f2 = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,x2,Absolute)
      f3 = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,x3,Absolute)
      Do While ((Abs(x4-x1) .gt. tol_x*(Abs(x1)+Abs(x2))) .and.
     &          (Abs(f3-f2) .gt. tol_f*(Abs(f2)+Abs(f3))))
        If (f2 .lt. f3) Then
          x4 = x3
          x3 = x2
          f3 = f2
          x2 = RM*x3 + Ratio*x1
          f2 = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,x2,Absolute)
        Else
          x1 = x2
          x2 = x3
          f2 = f3
          x3 = RM*x2 + Ratio*x4
          f3 = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,x3,Absolute)
        End If
      End Do
      If (f2 .lt. f3) Then
        xmin = x2
        Golden2 = f2
      Else
        xmin = x3
        Golden2 = f3
      End If

      End Function Golden2

#endif
