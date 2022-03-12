!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
!-- Rotate multipole.
!
      Subroutine Rotation_qmstat(iL,dMul,Rotte,Sigge)
      Implicit Real*8 (a-h,o-z)

      Parameter (MxMltp=2)

      Dimension dMul((MxMltp+1)*(MxMltp+2)/2),Rotte(3,3)
      Dimension dMTrans(6),TD(6,6)
#include "warnings.h"

!
!-- Charge, trivial to rotate.
!
      If(iL.eq.0) then
        dMul(1)=dMul(1)
!
!-- Dipole, transforms as a vector. Sigge controls that if the
!   multipole is located not in origin, but at the other end,
!   i.e. molecule A, then any odd occurrence of z should be
!   mirrored. Applies for the quadrupole as well, see below.
!
      ElseIf(iL.eq.1) then
        d1=dMul(1)
        d2=dMul(2)
        d3=dMul(3)
        dMul(1)=Rotte(1,1)*d1+Rotte(1,2)*d2+Rotte(1,3)*d3
        dMul(2)=Rotte(2,1)*d1+Rotte(2,2)*d2+Rotte(2,3)*d3
        dMul(3)=Rotte(3,1)*d1+Rotte(3,2)*d2+Rotte(3,3)*d3
        dMul(1)=dMul(1)
        dMul(2)=dMul(2)
        dMul(3)=Sigge*dMul(3)
!
!-- Quadrupole, transforms as a quadratic form. Also, transform
!   to spherical representation.
!
      ElseIf(iL.eq.2) then
!
!---- Compute the transformation matrix for second-moments.
!
        Call M2Trans(Rotte,TD)
!
!---- Transform. Sigge is explained above.
!
        Do i=1,6
          dMTrans(i)=0.0d0
          Do j=1,6
            dMTrans(i)=dMTrans(i)+TD(i,j)*dMul(j)
          Enddo
        Enddo
        Do i=1,6
          Sig=1.0d0
          If(i.eq.3.or.i.eq.5)Sig=Sigge
          dMul(i)=dMTrans(i)*Sig
        Enddo
!
!---- Go to spherical representation.
!
        Call Spherical(dMul)
      Else
        Write(6,*)'Nope!, Error in sl_grad'
        Call Quit(_RC_IO_ERROR_READ_)
      Endif

      Return
      End
