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
* Copyright (C) 1995, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine ShiftHess(Hess,shift,nDim,nDim2)
C!
C!  Purpose:
C!    Shifts Hessian to make it positive definite.
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
      Real*8 Hess (nDim,nDim2)
c       Real*8 U(nDim,nDim2)
c       Real*8 Hess_lowT( nDim*(nDim+1)/2 )
      Real*8  epsilon
      Real*8  eigen_min
      Logical  shift
#include "WrkSpc.fh"
      Call GetMem('U','Allo','Real',ipU,nDim*nDim2)
      Call GetMem('Hess_low','Allo','Real',
     &  ipHess_lowT,nDim*(nDim+1)/2)

C!
C!---- Initialize.
      nDim1 = nDim+1
      nDimSqr = nDim**2
C!
      k = 0
      Do i = 1,nDim
      Do j = 1,i
      k = k+1
      Work(ipHess_lowT+k-1) = Hess(i,j)
      End Do
      End Do
      call dcopy_(nDimSqr,0.0d0,0,Work(ipU),1)
      call dcopy_(nDim,1.0d0,0,Work(ipU),nDim1)
      Call Jacob(Work(ipHess_lowT),Work(ipU),nDim,nDim)
      Call Jacord(Work(ipHess_lowT),Work(ipU),nDim,nDim)
      eigen_min = Work(ipHess_lowT)
      shift = eigen_min.lt.0.0d0
      If ( shift ) Then
      epsilon = 2*eigen_min
      Do i = 1,nDim
      Hess(i,i) = Hess(i,i)-epsilon
      End Do
      End If
      Call GetMem('U','Free','Real',ipU,nDim*nDim2)
      Call GetMem('Hess_low','Free','Real',
     &  ipHess_lowT,nDim*(nDim+1)/2)
C!
      End
