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
* Copyright (C) 2002, Roland Lindh                                     *
************************************************************************
      Subroutine Assemble_dVdB(NAInt,EFInt,nZeta,la,lb,A,B,C)
************************************************************************
*                                                                      *
*     Object: to assemble the derivative of the nuclear attractoion    *
*             integrals with respect to the magnetic field.            *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics,                 *
*             University of Lund, SWEDEN                               *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 NAInt(nZeta*((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)),
     &       EFInt(nZeta*((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2),3),
     &       A(3), B(3), C(3), RAB(3)
*
      RAB(1)=A(1)-B(1)
      RAB(2)=A(2)-B(2)
      RAB(3)=A(3)-B(3)
*
*---- Recombine in place!
*
      nVec=nZeta*((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)
      Do iVec = 1, nVec
         EFInt_x=EFInt(iVec,1)
         EFInt_y=EFInt(iVec,2)
         EFInt_z=EFInt(iVec,3)
         EFInt(iVec,1)=RAB(2)*(EFInt_z+C(3)*NAInt(iVec))
     &                -RAB(3)*(EFInt_y+C(2)*NAInt(iVec))
         EFInt(iVec,2)=RAB(3)*(EFInt_x+C(1)*NAInt(iVec))
     &                -RAB(1)*(EFInt_z+C(3)*NAInt(iVec))
         EFInt(iVec,3)=RAB(1)*(EFInt_y+C(2)*NAInt(iVec))
     &                -RAB(2)*(EFInt_x+C(1)*NAInt(iVec))
      End Do
*
      Return
      End
