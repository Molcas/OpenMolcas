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
* Copyright (C) 1994, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine Freq(Hess,G,V,Lambda,B,Bnew,qMat,nOsc,NumOfAt)
C!
C!  Purpose:
C!    Find eigenvalues and eigenvectors of FG matrix.
C!    The eigenvalues are stored in the array Lambda and the eigen-
C!    vectors are stored as the columns of the matrix V.
C!    The eigenvectors are then used to calculate the cartesian
C!    displacements, which are stored in matrix X.
C!
C!  Input:
C!    Hess     : Real*8 two dimensional array -  contains
C!               the force constants expressed in internal
C!    G        : Real*8 two dimensional array.
C!    B        : Real*8 two dimensional array.
C!    Bnew     : Real*8 two dimensional array.
C!
C!  Output:
C!    V        : Real*8 two dimensional array  - contains
C!               the eigenvectors of F*G as columns.
C!    Lambda   : Real*8 array  - contains the eigenvalues
C!               of F*G.
C!    qMat     : Real*8 array - contains the cartesian
C!               displacements of the atoms.
C!
C!  Local:
C!    U,Tmp    : Real*8 two dimensional arrays.
C!
C!  Uses:
C!    LinAlg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1994.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 Hess(nOsc,nOsc),G(nOsc,nOsc),V(nOsc,nOsc)
      Real*8 B(3*NumOfAt,nOsc),Bnew(3*NumOfAt,nOsc),
     &  qMat(3*NumOfAt,nOsc)
      Real*8 Lambda(nOsc)
#include "WrkSpc.fh"

C!
C!---- Get dimensions.
      NumInt  = nOsc
C!
C!---- Solve secular equation.
      Call SolveSecEq(Hess,NumInt,V,G,Lambda)
C!
C!---- get memory space for U.
      Call GetMem('U','Allo','Real',ipU,nOsc*nOsc)
C!
C!---- Copy matrix containing eigenvectors to U, because this matrix
C!     will be destroyed when subroutine Dool is called.
      nSqrInt = NumInt**2
      call dcopy_(nSqrInt,V,1,Work(ipU),1)
C!
C!     Calculate the cartesian diplacements, i.e. solve
C!
C!              qMat = B * ( B(T) * B )^(-1) * V.
C!
C!---- get memory for matrix Temp.
      Call GetMem('Temp','Allo','Real',ipTemp,nOsc*nOsc)
C!
      Call DGEMM_('T','N',
     &            NumInt,NumInt,3*NumOfAt,
     &            1.0d0,B,3*NumOfAt,
     &            B,3*NumOfAt,
     &            0.0d0,Work(ipTemp),NumInt)
      Call Dool(Work(ipTemp),NumInt,NumInt,Work(ipU),
     &  NumInt,NumInt,Det)
      Call DGEMM_('N','N',
     &            3*NumOfAt,NumInt,NumInt,
     &            1.0d0,B,3*NumOfAt,
     &            Work(ipU),NumInt,
     &            0.0d0,qMat,3*NumOfAt)
C!
C!---- free memory space of Temp and U.
      Call GetMem('U','Free','Real',ipU,nOsc*nOsc)
      Call GetMem('Temp','Free','Real',ipTemp,nOsc*nOsc)
C!
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Bnew)
      End
