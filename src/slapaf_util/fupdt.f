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
      Subroutine FUpdt(nInter,FEq,gi_2,gi_1,g_ref,
     &                            qi_2,qi_1,q_ref,u,v,w)
************************************************************************
*                                                                      *
* Object: to do a rank-three update on the anharmonic constants.       *
*                                                                      *
*         Observe that SlapAf stores the forces rather than the        *
*         gradients.                                                   *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 FEq(nInter,nInter,nInter),
     &       gi_2(nInter), gi_1(nInter), g_ref(nInter),
     &       qi_2(nInter), qi_1(nInter), q_ref(nInter),
     &       u(nInter), v(nInter), w(nInter), lambda
*
*     Call RecPrt(' FEq',' ',FEq,nInter**2,nInter)
*     Call RecPrt('  gi-2',' ', gi_2,1,nInter)
*     Call RecPrt('  gi-1',' ', gi_1,1,nInter)
*     Call RecPrt('  g_ref',' ', g_ref,1,nInter)
*     Call RecPrt('  qi-2',' ', qi_2,1,nInter)
*     Call RecPrt('  qi-1',' ', qi_1,1,nInter)
*     Call RecPrt('  q_ref',' ', q_ref,1,nInter)
*
      Do i = 1, nInter
         u(i) = -(gi_1(i)-g_ref(i))
         v(i) = -(gi_2(i)-g_ref(i))
      End Do
      rLHS=DDot_(nInter,qi_2,1,u,1)-DDot_(nInter,q_ref,1,u,1)
     &    -DDot_(nInter,qi_1,1,v,1)+DDot_(nInter,q_ref,1,v,1)
      Write (6,*) 'FUpdt: LHS=',rLHS
      RHS=Zero
      Do i = 1, nInter
         Do k = 1, nInter
            Do l = 1, nInter
               RHS = RHS + FEq(i,k,l)*
     &               (qi_1(i)-q_ref(i))*
     &               (qi_2(k)-q_ref(k))*
     &               (qi_2(l)-qi_1(l))
            End Do
         End Do
      End Do
      RHS = Half*RHS
      Write (6,*) 'FUpdt: RHS=',RHS
      lambda=rLHS - RHS
      Write (6,*) ' FUpdt: lambda=',lambda
      Do i = 1, nInter
         w(i) = v(i) - u(i)
      End Do
      Call RecPrt('u',' ',u,1,nInter)
      Call RecPrt('v',' ',v,1,nInter)
      Call RecPrt('w',' ',w,1,nInter)
*
      ux = DDot_(nInter,u,1,qi_1,1)-DDot_(nInter,u,1,q_ref,1)
      uy = DDot_(nInter,u,1,qi_2,1)-DDot_(nInter,u,1,q_ref,1)
      vx = DDot_(nInter,v,1,qi_1,1)-DDot_(nInter,v,1,q_ref,1)
      vy = DDot_(nInter,v,1,qi_2,1)-DDot_(nInter,v,1,q_ref,1)
      wx = DDot_(nInter,w,1,qi_1,1)-DDot_(nInter,w,1,q_ref,1)
      wy = DDot_(nInter,w,1,qi_2,1)-DDot_(nInter,w,1,q_ref,1)
      lambda=Two*lambda/(ux*vy*(wy-wx)
     &                  +vx*wy*(uy-ux)
     &                  +wx*uy*(vy-vx))
      Write (6,*) ' FUpdt: lambda=',lambda
      Do i = 1, nInter
         Do k = 1, nInter
            Do l = 1, nInter
               FEq(i,k,l) = FEq(i,k,l) + lambda*(
     &                         u(i)*v(k)*w(l)
     &                       + v(i)*w(k)*u(l)
     &                       + w(i)*u(k)*v(l) )
            End Do
         End Do
      End Do
*
*
*---- Check on the third order condition
*
      Do i = 1, nInter
         u(i) = -(gi_1(i)-g_ref(i))
         v(i) = -(gi_2(i)-g_ref(i))
      End Do
      rLHS=DDot_(nInter,qi_2,1,u,1)-DDot_(nInter,q_ref,1,u,1)
     &    -DDot_(nInter,qi_1,1,v,1)+DDot_(nInter,q_ref,1,v,1)
      Write (6,*) 'FUpdt: LHS(qNR)=',rLHS
      RHS=Zero
      Do i = 1, nInter
         Do k = 1, nInter
            Do l = 1, nInter
               RHS = RHS + FEq(i,k,l)*
     &               (qi_1(i)-q_ref(i))*
     &               (qi_2(k)-q_ref(k))*
     &               (qi_2(l)-qi_1(l))
            End Do
         End Do
      End Do
      RHS = Half*RHS
      Write (6,*) 'FUpdt: RHS(qNR)=',RHS
*
      Return
      End
