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
*
* Compute and store roots and weights for a shifted Legendre quadrature,
* used for boot-strapping Rys roots and weights.
* Different sets of roots and weights are computed
*
      Module Leg_RW
      Implicit None
      Integer, Dimension(11), Parameter :: naux=[30,35,40,45,50,
     &                                           55,60,65,70,75,300]
      Real*8, Dimension(:,:), Allocatable :: Leg_r, Leg_w

      Contains

      Subroutine SetAux(eps)
      Real*8, Intent(In) :: eps
      Integer, Parameter :: nquad=Size(naux)
      Real*8, Dimension(:), Allocatable :: a, b
      Integer :: maux, i, j, Err
#include "stdalloc.fh"
#include "real.fh"
      If (Allocated(Leg_r)) Return
      maux = MaxVal(naux)
      Call mma_allocate(Leg_r,maux,nquad,label="Leg_r")
      Call mma_allocate(Leg_w,maux,nquad,label="Leg_w")
      Call mma_allocate(a,maux)
      Call mma_allocate(b,maux)
      Do j = 1, nquad
        Do i = 1, naux(j)
          a(i) = Half
          If (i == 1) Then
            b(1) = One
          Else
            b(i) = Quart/(Four-One/(i-1)**2)
          End If
        End Do
        Call GaussQuad(naux(j),a,b,eps,Leg_r(1,j),Leg_w(1,j),Err)
        If (Err.ne.0) Then
          write(6,*) Err
          Call WarningMessage(2,'Error in GaussQuad')
          Call AbEnd()
        End If
        Do i = 1, naux(j)
          Leg_r(i,j)=Leg_r(i,j)*Leg_r(i,j)
        End Do
      End Do
      Call mma_deallocate(a)
      Call mma_deallocate(b)
      End Subroutine SetAux

      Subroutine UnSetAux
#include "stdalloc.fh"
      If (Allocated(Leg_r)) Call mma_deallocate(Leg_r)
      If (Allocated(Leg_w)) Call mma_deallocate(Leg_w)
      End Subroutine UnSetAux

      End Module Leg_RW
