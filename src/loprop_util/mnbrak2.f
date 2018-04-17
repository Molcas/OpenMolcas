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

      Subroutine MnBrak2(ax,bx,cx,fa,fb,fc,f,
     &                 q_a,q_b,dipole_a,dipole_b,r_a,r_b)
      Implicit None
      Real*8 :: ax, bx, cx, fa, fb, fc, vx, fv, coefA, coefB
      Logical :: Def
      Real*8, Parameter :: Ratio = 0.5D0*(Sqrt(5.0D0)+1.0D0)
      Real*8, Parameter :: Thr = 1.0D-20, Lim = 100.0D0
c External function f and its arguments
      Real*8, External :: f
      Real*8 :: q_a,q_b,dipole_a,dipole_b,r_a,r_b
      Logical, Parameter :: Absolute=.True.

#include "real.fh"
      fa = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,ax,Absolute)
      fb = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,bx,Absolute)
      If (fa .lt. fb) Then
        cx = ax
        ax = bx
        bx = cx
        fc = fa
        fa = fb
        fb = fc
      End If
      ! three points such that b is between a and c,
      ! and f(a) > f(b) > f(c), stop when f(c) > f(b)
      cx = bx + Ratio*(bx-ax)
      fc = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,cx,Absolute)
      Do While (fc .le. fb)
        write(6,*) ax,bx,cx
        Def = .True.
        ! try a parabolic fitting
        coefA = (   cx*(fb-fa)+   bx*(fa-fc)+   ax*(fc-fb))
        coefB = (cx**2*(fa-fb)+bx**2*(fc-fa)+ax**2*(fb-fc))
        ! only worth it if the points are not linear and
        ! if the 2nd derivative is positive
        If ((Abs(coefA) .gt. Thr) .and. coefA*(ax-cx) .gt. Zero) Then
          vx = -Half*coefB/coefA
          ! v is between b and c
          If ((cx-vx)*(vx-bx) .gt. Zero) Then
            fv = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,vx,Absolute)
            ! minimum between b and c
            If (fv .lt. fc) Then
              ax = bx
              bx = vx
              fa = fb
              fb = fv
              Return
            ! minimum between a and v
            Else If (fv .gt. fb) Then
              cx = vx
              fc = fv
              Return
            End If
          ! v is beyond c, but within limits
          Else If ((bx+Lim*(cx-bx)-vx)*(vx-cx) .gt. Zero) Then
            fv = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,vx,Absolute)
            ! whatever happens, replace c,v -> b,c
            bx = cx
            cx = vx
            fb = fc
            fc = fv
            ! minimum between b and v
            If (fv .gt. fc) Then
              ax = bx
              fa = fb
              Return
            End If
          ! v is beyond the limit, cutoff to the limit
          Else If ((vx-cx)*(cx-bx) .gt. Zero) Then
            fv = bx + Lim*(cx-bx)
            Def = .False.
          End If
        End If
        ! unless the fit went beyond limits, use default step
        If (Def) Then
          vx = cx + Ratio*(cx-bx)
          fv = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,vx,Absolute)
        End If
        ax = bx
        bx = cx
        cx = vx
        fa = fb
        fb = fc
        fc = fv
      End Do

      End Subroutine MnBrak2
