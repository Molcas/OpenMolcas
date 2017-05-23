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
* Copyright (C) 2000, Gunnar Karlstrom                                 *
*               2000, Roland Lindh                                     *
************************************************************************
      Subroutine qlm(gx,gy,gz,qa,dax,day,daz,lmax_,Cavxyz)
************************************************************************
*                                                                      *
*     Object: to reexpand the charge and the dipole moment at a given  *
*             point as a multipole moment expansion at origin.         *
*                                                                      *
*     Authors: G. Karlstroem                                           *
*              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
*                                                                      *
*              and                                                     *
*                                                                      *
*              R. Lindh                                                *
*              Dept. of Chem. Phys., Univ. of Lund, Sweden.            *
*                                                                      *
*              March 2000                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Cavxyz((lMax_+1)*(lMax_+2)*(lMax_+3)/6)
*
*---- Statement function
*
      iOff(ix,iy,iz)=(ix+iy+iz)*(ix+iy+iz+1)*(ix+iy+iz+2)/6
      Index(ix,iy,iz) = iOff(ix,iy,iz) + (iy+iz)*(iy+iz+1)/2 + iz + 1
*
      lmax = lmax_ - 1
*
      Do ix=0,lmax
         If (ix.eq.0) Then
            xeff = One
         Else
            xeff = gx**ix
         End If
         ax=DBLE(ix)+One
         Do iy=0,lmax-ix
            If (iy.eq.0) Then
               xyeff = xeff
            Else
               xyeff = xeff* gy**iy
            End If
            ay=DBLE(iy)+One
            Do iz=0,lmax-ix-iy
               If (iz.eq.0) Then
                  xyzeff = xyeff * gz**iz
               Else
                  xyzeff = xyeff * gz**iz
               End If
               az=DBLE(iz)+One
*
*------------- Charge term
*
               Cavxyz(Index(ix,iy,iz)) = xyzeff*qa
     &                                 + Cavxyz(Index(ix,iy,iz))
*
*------------- Dipole terms
*
               Cavxyz(Index(ix+1,iy,iz)) = xyzeff*dax*ax
     &                                   + Cavxyz(Index(ix+1,iy,iz))
               Cavxyz(Index(ix,iy+1,iz)) = xyzeff*day*ay
     &                                   + Cavxyz(Index(ix,iy+1,iz))
               Cavxyz(Index(ix,iy,iz+1)) = xyzeff*daz*az
     &                                   + Cavxyz(Index(ix,iy,iz+1))
*
            End Do
         End Do
      End Do
*
      Return
      End
