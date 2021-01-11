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
*  Branching_Plane_Update
*
*> @brief
*>   Apply the Branching plane update method of Maeda et al.
*> @author Roland Lindh
*> @modified_by Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> The branching plane of conical intersection is defined by two vectors,
*> the gradient difference vector (GDV) and the coupling derivative vector
*> (CDV). Equivalently, it can be defined with two orthonormal vectors
*> \f$ x \f$ and \f$ y \f$, such that \f$ x \f$ has the same direction as
*> GDV and \f$ \{x,y\} \f$ spans the same space as GDV and CDV.
*>
*> This routine obtains an approximate \f$ y \f$ vector, by assuming that
*> the branching plane defined by \f$ \{x,y\} \f$ changes the least from
*> iteration to iteration \cite Mae2010-JCTC-6-1538. Thus:
*>
*> \f[ y_k = \beta y_{k-1} - \alpha x_{k-1} \\
*>     \alpha = \frac{y_{k-1}\cdot x_k}{w} \\
*>     \beta  = \frac{x_{k-1}\cdot x_k}{w} \\
*>     w      = \sqrt{(y_{k-1}\cdot x_k)^2+(x_{k-1}\cdot x_k)^2} \f]
*>
*> As an initial guess for the CDV, the average gradient vector (AGV) at
*> the first iteration is used.
*>
*> @param[in]     AGV   Average gradient vector(s)
*> @param[in]     DGV   Difference gradient vector(s)
*> @param[in,out] CDV   Approximate coupling derivative vector, \f$ y \f$
*> @param[in]     n     Size of the vectors
*> @param[in]     nIter Iteration number
************************************************************************
      Subroutine Branching_Plane_Update(AGV,DGV,CDV,n,nIter)
      Implicit Real*8 (a-h,o-z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8 AGV(n,nIter), DGV(n,nIter), CDV(n)
      Real*8, Allocatable:: x0(:), x1(:)
*
      iRout = 31
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.6) Then
         Write (6,*) 'Branching plane'
         Write (6,*) 'n,nIter=',n,nIter
         Call RecPrt('AGV',' ',AGV,n,nIter)
         Call RecPrt('DGV',' ',DGV,n,nIter)
         Call RecPrt('CDV (init)',' ',CDV,n,1)
      End If
      Call mma_allocate(x0,n,Label='x0')
      Call mma_allocate(x1,n,Label='x1')
*
*     Get the directional vector for the first difference gradient
*     vector (DGV).
*
      call dcopy_(n,DGV,1,x0,1)
      r = One / Sqrt(DDot_(n,x0,1,x0,1))
      Call DScal_(n,r,x0,1)
      call dcopy_(n,x0,1,x1,1)
*
*     Initial guess for the coupling derivative vector (CDV) is the
*     average gradient vector (AGV).
*     ... or rather its component perpendicular to DGV
*     so, remove the projection of CDV along DGV and renormalize
*
      call dcopy_(n,AGV,1,CDV,1)
      proj = DDot_(n,CDV,1,x0,1)
      Call DaXpY_(n,-proj,x0,1,CDV,1)
      r = One / Sqrt(DDot_(n,CDV,1,CDV,1))
      Call DScal_(n,r,CDV,1)
      If (iPrint.ge.6) Then
         Call RecPrt('CDV(0)',' ',CDV,n,1)
      End If
*
*     Apply the MOM update method to correct the guessed CDV.
*
*     ipx0: xk-1, ~DGV of previous iteration
*     ipx1: xk,   ~DGV of current iteration
*
      Do iter=2,nIter
         call dcopy_(n,DGV(1,iter),1,x1,1)
         r = One / Sqrt(DDot_(n,x1,1,x1,1))
         Call DScal_(n,r,x1,1)
*
         xx=DDot_(n,x0,1,x1,1)
         yx=DDot_(n,CDV,1,x1,1)
         yx_xx=Sqrt(yx**2+xx**2)
*
*        different signs from the paper, should not matter,
*        but will keep the y vector sign if x does not change
         alpha=-yx/yx_xx
         beta=xx/yx_xx
         Call DScal_(n,beta,CDV,1)
         Call DaXpY_(n,alpha,x0,1,CDV,1)
*
*        remove the projection of CDV along DGV and renormalize
*
         If (iPrint.ge.6) Then
            Write (6,*)
            Write (6,*) 'iter=',iter
            Write (6,*) 'r(DGV)=',r
            Write (6,*) 'xx=',xx
            Write (6,*) 'yx=',yx
            Write (6,*) 'alpha,beta=',alpha,beta
         End If
*
         proj=DDot_(n,CDV,1,x1,1)
         Call DaXpY_(n,-proj,x1,1,CDV,1)
         r = One / Sqrt(DDot_(n,CDV,1,CDV,1))
         Call DScal_(n,r,CDV,1)
*
         If (iPrint.ge.6) Then
            Write (6,*) 'r(CDV)=',r
         End IF
*
         If (iter.ne.nIter) call dcopy_(n,x0,1,x1,1)
*
      End Do
*
      Call mma_deallocate(x1)
      Call mma_deallocate(x0)
      If (iPrint.ge.6) Then
         Call RecPrt('CDV',' ',CDV,n,1)
      End If
*
      Return
      End
