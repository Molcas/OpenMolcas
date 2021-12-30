!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************
      Subroutine Rot_st(cMO_s,cMO_t,nBasis,Gamma_rot,Debug)
!
!     Author: Y. Carissan.
!
      Implicit Real*8(a-h,o-z)
#include "real.fh"
      Real*8 cMO_s(nBasis),cMO_t(nBasis)
      Logical Debug
!
      If (Gamma_rot.eq.Zero) Return
!
      cosGamma_rot=cos(Gamma_rot)
      sinGamma_rot=sin(Gamma_rot)
      If (Debug) Then
        Write(6,*) 'cos(Gamma)=',cosGamma_rot
        Write(6,*) 'sin(Gamma)=',sinGamma_rot
      End If
!
      Do iBas=1,nBasis
        cs=cMO_s(iBas)
        ct=cMO_t(iBas)
        cMO_s(iBas)= cosGamma_rot*cs + sinGamma_rot*ct
        cMO_t(iBas)=-sinGamma_rot*cs + cosGamma_rot*ct
      End Do
!
      Return
      End
