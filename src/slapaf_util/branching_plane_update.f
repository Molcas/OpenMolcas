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
      Subroutine Branching_Plane_Update(AGV,DGV,CDV,n,nIter)
************************************************************************
*     Apply the Braching plane update method of S. Maeda, K. Ohno and  *
*     K. Morokuma, JCTC, 6 (2010), 1538, DOI:10.1021/ct1000268         *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 AGV(n,nIter), DGV(n,nIter), CDV(n)
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
      Call Allocate_Work(ipx0,n)
      Call Allocate_Work(ipx1,n)
*
*     Get the directional vector for the first difference gradient
*     vector (DGV).
*
      call dcopy_(n,DGV,1,Work(ipx0),1)
      r = One / Sqrt(DDot_(n,Work(ipx0),1,Work(ipx0),1))
      Call DScal_(n,r,Work(ipx0),1)
      call dcopy_(n,Work(ipx0),1,Work(ipx1),1)
*
*     Initial guess for the coupling derivative vector (CDV) is the
*     average gradient vector (AGV).
*     ... or rather its component perpendicular to DGV
*     so, remove the projection of CDV along DGV and renormalize
*
      call dcopy_(n,AGV,1,CDV,1)
      proj = DDot_(n,CDV,1,Work(ipx0),1)
      Call DaXpY_(n,-proj,Work(ipx0),1,CDV,1)
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
         call dcopy_(n,DGV(1,iter),1,Work(ipx1),1)
         r = One / Sqrt(DDot_(n,Work(ipx1),1,Work(ipx1),1))
         Call DScal_(n,r,Work(ipx1),1)
*
         xx=DDot_(n,Work(ipx0),1,Work(ipx1),1)
         yx=DDot_(n,CDV,1,Work(ipx1),1)
         yx_xx=Sqrt(yx**2+xx**2)
*
*        different signs from the paper, should not matter,
*        but will keep the y vector sign if x does not change
         alpha=-yx/yx_xx
         beta=xx/yx_xx
         Call DScal_(n,beta,CDV,1)
         Call DaXpY_(n,alpha,Work(ipx0),1,CDV,1)
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
         proj=DDot_(n,CDV,1,Work(ipx1),1)
         Call DaXpY_(n,-proj,Work(ipx1),1,CDV,1)
         r = One / Sqrt(DDot_(n,CDV,1,CDV,1))
         Call DScal_(n,r,CDV,1)
*
         If (iPrint.ge.6) Then
            Write (6,*) 'r(CDV)=',r
         End IF
*
         If (iter.ne.nIter) call dcopy_(n,Work(ipx0),1,Work(ipx1),1)
*
      End Do
*
      Call Free_Work(ipx1)
      Call Free_Work(ipx0)
      If (iPrint.ge.6) Then
         Call RecPrt('CDV',' ',CDV,n,1)
      End If
*
      Return
      End
