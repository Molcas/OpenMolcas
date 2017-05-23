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
      Subroutine hmod(gx,gy,gz,V,EFx,EFy,EFz,Cavxyz,lmax_)
************************************************************************
*                                                                      *
*     Object: to compute the potential and electric field in a point   *
*             due to the multipole moment expansion at origin.         *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Cavxyz((lMax_+1)*(lMax_+2)*(lMax_+3)/6)
*
*---- Statement function
*
      iOff(ix,iy,iz)=(ix+iy+iz)*(ix+iy+iz+1)*(ix+iy+iz+2)/6
      Index(ix,iy,iz) = iOff(ix,iy,iz) + (iy+iz)*(iy+iz+1)/2 + iz + 1
*
      V  = Zero
      EFx= Zero
      EFy= Zero
      EFz= Zero
*
*     nCavxyz=(lMax_+1)*(lMax_+2)*(lMax_+3)/6
*     Write (*,*) 'lMax_,nCavxyz=',lMax_,nCavxyz
      lmax=lmax_-1
      Do ix = 0, lmax
         If (ix.eq.0) Then
            xeff = One
         Else
            xeff = gx**ix
         End If
         ax = DBLE(ix) + One
         Do iy = 0, lmax-ix
            If (iy.eq.0) Then
               xyeff = xeff
            Else
               xyeff = xeff * gy**iy
            End If
            ay = DBLE(iy) + One
            Do iz = 0, lmax-ix-iy
               If (iz.eq.0) Then
                  xyzeff = xyeff
               Else
                  xyzeff = xyeff * gz**iz
               End If
               az = DBLE(iz) + One
*              Write (*,*) ix,iy,iz,Index(ix,iy,iz),
*    &                     Index(ix+1,iy,iz),Index(ix,iy+1,iz),
*    &                     Index(ix,iy,iz+1)
*
*------------- Charge term
*
               V=V+xyzeff*Cavxyz(Index(ix  ,iy  ,iz  ))
*
*------------- Dipole terms
*
               EFx=EFx+ax*xyzeff*Cavxyz(Index(ix+1,iy  ,iz  ))
               EFy=EFy+ay*xyzeff*Cavxyz(Index(ix  ,iy+1,iz  ))
               EFz=EFz+az*xyzeff*Cavxyz(Index(ix  ,iy  ,iz+1))
*
            End Do
         End Do
      End Do
*
      Return
      End
