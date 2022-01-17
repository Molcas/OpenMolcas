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
      Subroutine Do_GGL(L_Eff,mPt,ipR)
************************************************************************
*                                                                      *
*     Computes datas useful for the angular quadrature.                *
*                                                                      *
************************************************************************
      use nq_Grid, only: Pax
      Implicit Real*8 (a-h,o-z)
#include "nq_info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: Th(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate angular grid from Gauss and Gauss-Legendre quadrature
*
*---- Theta (polar angle): 0 =< theta =< pi
*     Gauss-Legendre Quadrature (L_Quad+1)/2 points
*---- Phi (azimuthal angle): 0=< phi =< 2*pi
*     Gauss-Quadrature (L_Quad+1) points
*
      nTheta = (L_Eff+1)/2
      nPhi   = L_Eff+1
*
      mPt    = nTheta*nPhi
      Call GetMem('AngRW','Allo','Real',ipR,4*mPt)
*
      Call mma_allocate(Th,2,nTheta,Label='Th')
*
      Call GauLeg(-One,One,Th,nTheta)
*
      iOff = ipR
      Do iTheta = 1, nTheta
         Cos_Theta = Th(1,iTheta)
         w_Theta   = Th(2,iTheta)
         Sin_Theta = Sqrt(One-Cos_Theta**2)
*
         Do iPhi = 1, nPhi
            Call Phi_point(iPhi,nPhi,Cos_Phi,Sin_Phi,w_Phi)
*
            x = Sin_Theta*Cos_Phi
            y = Sin_Theta*Sin_Phi
            z = Cos_Theta
            Work(iOff  )=Pax(1,1)*x+Pax(1,2)*y+Pax(1,3)*z
            Work(iOff+1)=Pax(2,1)*x+Pax(2,2)*y+Pax(2,3)*z
            Work(iOff+2)=Pax(3,1)*x+Pax(3,2)*y+Pax(3,3)*z
            Work(iOff+3)=w_Theta*w_Phi
            iOff = iOff + 4
*
         End Do    ! iPhi
*
      End Do       ! iTheta
*
      Call mma_deallocate(Th)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
