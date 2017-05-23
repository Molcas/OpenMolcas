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
* Copyright (C) 2005, Christian Ander                                  *
************************************************************************
      Subroutine ts_bfgs(dq,y,gi,H,nH)
      Implicit None
*     Hessian update method; TS-BFGS ; Bofill - "Remarks on the
*     Updated Hessian Matrix Methods" 2003.
*
*     Vectors in the article are column vectors.
*
*     Implemented by Christian Ander, 2005, christian@eximius.se
*                                                                      *
************************************************************************
*                                                                      *
*     gi    :  gradient                   (nH)
*     dq    :  Perturbation Geometry      (nH)
*     y     :  Difference in Gradient     (nH)
*     WorkM :  Temporary working matrix   (nH,nH)
*     WorkV :  Temporary working vector   (nH)
*     WorkR :  Temporary working variable (real*8)
*     H     :  Hessian                    (nH,nH)
*     nH    :  Hessian size               (integer)
*     f,a,b :  Multi-used variables       (real*8)
*     v,u   :  Multi-used vectors         (nH)
*
#include "real.fh"
#include "stdalloc.fh"
      Integer nH, i, j
      Real*8 H(nH,nH), dq(nH), y(nH), gi(nH)
      Real*8 a, b, f, WorkR, ddot_
      Real*8, Dimension(:,:), Allocatable :: WorkM
      Real*8, Dimension(:), Allocatable :: WorkV, v, u
*
      Call mma_allocate(WorkM,nH,nH,label="WorkM")
      Call mma_allocate(WorkV,nH,label="WorkV")
      Call mma_allocate(v,nH,label="v")
      Call mma_allocate(u,nH,label="u")
*                                                                      *
************************************************************************
*                                                                      *
*     Equation 23 calculates the new Hessian, some simplifications :
*
*     H_k+1 = H_k + 1/f (v*u^T + u*v^T - (y^T*dq - dq^T*B*dq)/f * u*u^T)
*
*     Where f = (y^T*dq)^2 + (dq^T|B|dq)^2 = a^2 + b^2
*     ,     u = y^T*dq*y + (dq^T|B|dq)|B|dq = a*y + b|B|dq
*     and   v = y - B*dq (quasi-Newton condition)
*
*define _DEBUG_
#ifdef _DEBUG_
*     Make a comment in logfile
      write(6,*) 'hello from ts_bfgs.f'
#endif
*
*---- Calculation of u = y^T*dq*y + (dq^T|B|dq)|B|dq
*
*     a = y^T * dq
*
      a = DDot_(nH,y,1,dq,1)
*
*     u = y^T*dq*y = a * y
*
      call dcopy_(nH,y,1,u,1)
      Call DScal_(nH,a,u,1)
*
*     WorkM = |B|
*
      Do j = 1,nH
         Do i = 1,nH
            WorkM(i,j) = Abs(H(i,j))
         End Do
      End Do
*
*     WorkV = |B|dq
*
      Call dGeMV_('N',nH,nH,
     &           One,WorkM,nH,
     &           dq,1,
     &           Zero,WorkV,1)
*
*     b = (dq^T|B|dq) = dq^T * WorkV
*
      b = DDot_(nH,dq,1,WorkV,1)
*
*     u = y^T*dq*y + (dq^T|B|dq)|B|dq = u + b * WorkV
*
      Call DaxPy_(nH,b,WorkV,1,u,1)
*
*---- Calculation of f = (y^T*dq)^2 + (dq^T|B|dq)^2 = a^2 + b^2
*
      f = a**2 + b**2
*
*---- Calculation of v = y - B * dq (quasi-Newton condition)
*
      call dcopy_(nH,y,1,v,1)
      Call dGeMV_('N',nH,nH,
     &           -One,H,nH,
     &           dq,1,
     &           One,v,1)
*
*---- Calculation of equation 23, the new Hessian
*
*     H_k+1 = H_k + 1/f (v*u^T + u*v^T - (a - dq^T*B*dq)/f * u*u^T)
*
*     WorkM = u*u^T
*
      Call DGEMM_('N','N',
     &            nH,nH,1,
     &            One,u,nH,
     &            u,1,
     &            Zero,WorkM,nH)
*
*     WorkV = dq^T*B
*
      Call DGEMM_('N','N',
     &            1,nH,nH,
     &            One,dq,1,
     &            H,nH,
     &            Zero,WorkV,1)
*
*     WorkR = (a - dq^T*B*dq)/f = (a - WorkV * dq) / f
*
      WorkR = (a - DDot_(nH,WorkV,1,dq,1)) / f
*
*     H_k+1 = H_k + 1/f ( v*u^T + u*v^T - (a - dq^T*B*dq)/f * u*u^T ) =
*             H_k + 1/f ( v*u^T + u*v^T - WorkR * WorkM )
*
      Do j = 1,nH
         Do i = 1,nH
            H(i,j) = H(i,j) + One/f * ( v(j)*u(i) + u(j)*v(i) -
     &               WorkR * WorkM(i,j) )
         End Do
      End Do
#ifdef _DEBUG_
*
*     Checking the quasi-Newton condition.
*
*     WorkV = H * dq
*
      Call dGeMV_('N',nH,nH,
     &           One,H,nH,
     &           dq,1,
     &           Zero,WorkV,1)
*
*     WorkV = WorkV - y = 0.00
*
      Call DaxPy_(nH,-One,y,1,WorkV,1)
      Call RecPrt('quasi-Newton',' ',WorkV,1,nH)
*
      write(6,*) 'good-bye from ts_bfgs.f'
#endif
*
      Call mma_deallocate(WorkM)
      Call mma_deallocate(WorkV)
      Call mma_deallocate(v)
      Call mma_deallocate(u)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(gi)
      End
