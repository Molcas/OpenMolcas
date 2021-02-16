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
      Subroutine qvar_to_var(var,x,grad,Hess,D3,D4,ref,qref,trfName,
     &  alpha, max_term,ndata,nvar)
C!
C!  Purpose:
C!    Tranform coordinates, gradient, Hessian, third derivatives and
C!    fourth derivates back to the coordinates originally specified.
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
      Implicit Real*8 ( a-h,o-z )
      Real*8 var (ndata,nvar)
      Real*8 x (nvar)
      Real*8 q (nvar)
      Real*8 t(nvar),u(nvar),v(nvar),s(nvar)
      Real*8 alpha (nvar)
      Real*8 grad(nvar),tempgrad(nvar)
      Real*8 Hess (nvar,nvar)
      Real*8 Temp (nvar,nvar)
      Real*8 D3  (nvar,nvar,nvar)
      Real*8 D3trans (nvar,nvar,nvar)
      Real*8 D4 (nvar,nvar,nvar,nvar)
      Real*8 D4trans (nvar,nvar,nvar,nvar)
      Real*8 ref(nvar),qref (nvar)
c       Integer nOrd( 4 )
      Character*80 trfName (nvar)
      Character*32  trfCode
      Logical  ijEq,jkEq,klEq
C!
C!
      do iv=1,nvar
      x(iv) = x(iv)+qref(iv)
      enddo
      Do ivar = 1,nvar
      trfcode = trfName(ivar)(1:32)
      ix = index(trfcode,'AS IT IS')
      ia = index(trfcode,'-AVG')
      ie = index(trfcode,'EXP')
      ir = index(trfcode,'RAD')
      id = index(trfcode,'DEG')
      ic = index(trfcode,'COS')
      is = index(trfcode,'SIN')
      If ( ic.gt.0 ) Then
      If ( ia.gt.0 ) Then
      x(ivar) = x(ivar)-ref(ivar)
      End If
      q(ivar) = acos(-x(ivar))
      t(ivar) = sin(q(ivar))
      u(ivar) = cos(q(ivar))
      v(ivar) =-sin(q(ivar))
      s(ivar) =-cos(q(ivar))
      Else If ( is.gt.0 ) Then
      If ( ia.gt.0 ) Then
      x(ivar) = x(ivar)-ref(ivar)
      End If
      q(ivar) = asin(-x(ivar))
      t(ivar) =-cos(q(ivar))
      u(ivar) = sin(q(ivar))
      v(ivar) = cos(q(ivar))
      s(ivar) =-sin(q(ivar))
      Else If ( ie.gt.0 ) Then
      q(ivar) = -log(1.0d0-x(ivar))/alpha(ivar)
      t(ivar) = alpha(ivar)*(1.0d0-x(ivar))
      u(ivar) =-alpha(ivar)**2*(1.0d0-x(ivar))
      v(ivar) = alpha(ivar)**3*(1.0d0-x(ivar))
      s(ivar) =-alpha(ivar)**4*(1.0d0-x(ivar))
      If ( ia.gt.0 ) Then
      q(ivar) = q(ivar)+ref(ivar)
      End If
      Else If ( ix.gt.0 ) Then
      If ( ia.gt.0 ) Then
      x(ivar) = x(ivar)+ref(ivar)
      End If
      q(ivar) = x(ivar)
      t(ivar) = 1.0d0
      u(ivar) = 0.0d0
      v(ivar) = 0.0d0
      s(ivar) = 0.0d0
      Else
      Write(6,*)' TRFCODE ERROR.'
      call abend()
      End if
      End Do
C!
C!---- Transform gradient.
      Do i = 1,nvar
      tempgrad(i) = t(i)*grad(i)
      End Do
C!
C!---- Transform Hessian.
      Do i = 1,nvar
      Do j = i,nvar
      If ( i.eq.j ) Then
      Temp(i,i) = t(i)**2*Hess(i,i)+u(i)*grad(i)
      Else
      Temp(i,j) = t(i)*t(j)*Hess(i,j)
      End If
      End Do
      End Do
      Do i = 1,nvar
      Do j = i,nvar
      Temp(j,i) = Temp(i,j)
      End Do
      End Do
C!
C!---- Transform third derivatives.
      Do i = 1,nvar
      Do j = i,nvar
      ijEq = ( i.eq.j )
      Do k = j,nvar
      jkEq = ( j.eq.k )
      If ( ijEq.and.jkEq ) Then
      D3trans(i,i,i) = t(i)*t(i)*t(i)*D3(i,i,i)+
     &                        3.0d0*u(i)*t(i)*Hess(i,i)+v(i)*grad(i)
      Else If ( ijEq.and.( .not.jkEq )) Then
      D3trans(i,i,k) = t(i)*t(i)*t(k)*D3(i,i,k)+
     &                        u(i)*t(k)*Hess(i,k)
      Else If (( .not.ijEq ).and.jkEq ) Then
      D3trans(i,j,j) = t(i)*t(j)*t(j)*D3(i,j,j)+
     &                        t(i)*u(j)*Hess(i,j)
      Else If (( .not.ijEq ).and.( .not.jkEq )) Then
      D3trans(i,j,k) = t(i)*t(j)*t(k)*D3(i,j,k)
      End If
      End Do
      End Do
      End Do
      Do i = 1,nvar
      Do j = i,nvar
      Do k = j,nvar
      D3trans(i,k,j) = D3trans(i,j,k)
      D3trans(j,i,k) = D3trans(i,j,k)
      D3trans(j,k,i) = D3trans(i,j,k)
      D3trans(k,i,j) = D3trans(i,j,k)
      D3trans(k,j,i) = D3trans(i,j,k)
      End Do
      End Do
      End Do
C!
C!---- Transform fourth derivatives.
      Do i = 1,nvar
      Do j = i,nvar
      ijEq = ( i.eq.j )
      Do k = j,nvar
      jkEq = ( j.eq.k )
      Do l = k,nvar
      klEq = ( k.eq.l )
      If ( ijEq.and.jkEq.and.klEq ) Then
      D4trans(i,i,i,i) = t(i)**4*D4(i,i,i,i)+
     &                           6.0d0*u(i)*t(i)**2*D3(i,i,i)+
     &                           4.0d0*v(i)*t(i)*Hess(i,i)+
     &                           3.0d0*u(i)*u(i)*Hess(i,i)+s(i)*grad(i)
      Else If ( ijEq.and.jkEq.and.( .not.klEq )) Then
      D4trans(i,i,i,l) = (t(i)*t(i)*t(i)*D4(i,i,i,l)+
     &                           3.0d0*u(i)*t(i)*D3(i,i,l)+v(i)*
     &                           Hess(i,l))*t(l)
      Else If (( .not.ijEq ).and.jkEq.and.klEq ) Then
      D4trans(i,j,j,j) = t(i)*(t(j)*t(j)*t(j)*
     &                D4(i,j,j,j)+
     &                3.0d0*u(j)*t(j)*D3(i,j,j)+v(j)*Hess(i,j))
      Else If (( .not.ijEq).and.jkEq.and.
     &             ( .not.klEq )) Then
      D4trans(i,j,j,l) = t(i)*(t(j)*t(j)*D4(i,j,j,l)+
     &                           u(j)*D3(i,j,l))*t(l)
      Else If ( ijEq.and.( .not.jkEq ).and.
     &             ( .not.klEq )) Then
      D4trans(i,i,k,l) = (t(i)*t(i)*D4(i,i,k,l)+
     &                           u(i)*D3(i,k,l))*t(k)*t(l)
      Else If (( .not.ijEq ).and.
     &             ( .not.jkEq ).and.klEq ) Then
      D4trans(i,j,k,k) = (t(k)*t(k)*D4(i,j,k,k)+
     &                           u(k)*D3(i,j,k))*t(i)*t(j)
      Else If ( ijEq.and.klEq ) Then
      D4trans(i,i,k,k) = t(i)*t(i)*t(k)*t(k)*
     &                D4(i,i,k,k)+
     &                t(i)*t(i)*u(k)*D3(i,i,k)+u(i)*t(k)*t(k)*
     &                D3(i,k,k)+
     &                u(i)*u(k)*Hess(i,k)
      Else If (( .not.ijEq ).and.( .not.jkEq ).and.
     &             ( .not.klEq )) Then
      D4trans(i,j,k,l) = t(i)*t(j)*t(k)*t(l)*
     &                D4(i,j,k,l)
      End If
      End Do
      End Do
      End Do
      End Do
      Do i = 1,nvar
      Do j = i,nvar
      Do k = j,nvar
      Do l = k,nvar
      D4trans(i,k,j,l) = D4trans(i,j,k,l)
      D4trans(i,j,l,k) = D4trans(i,j,k,l)
      D4trans(i,l,k,j) = D4trans(i,j,k,l)
      D4trans(i,l,j,k) = D4trans(i,j,k,l)
      D4trans(i,k,l,j) = D4trans(i,j,k,l)
C!
      D4trans(j,i,k,l) = D4trans(i,j,k,l)
      D4trans(j,k,i,l) = D4trans(i,j,k,l)
      D4trans(j,i,l,k) = D4trans(i,j,k,l)
      D4trans(j,l,k,i) = D4trans(i,j,k,l)
      D4trans(j,l,i,k) = D4trans(i,j,k,l)
      D4trans(j,k,l,i) = D4trans(i,j,k,l)
C!
      D4trans(k,i,j,l) = D4trans(i,j,k,l)
      D4trans(k,j,i,l) = D4trans(i,j,k,l)
      D4trans(k,i,l,j) = D4trans(i,j,k,l)
      D4trans(k,l,j,i) = D4trans(i,j,k,l)
      D4trans(k,l,i,j) = D4trans(i,j,k,l)
      D4trans(k,j,l,i) = D4trans(i,j,k,l)
C!
      D4trans(l,i,j,k) = D4trans(i,j,k,l)
      D4trans(l,j,i,k) = D4trans(i,j,k,l)
      D4trans(l,k,j,i) = D4trans(i,j,k,l)
      D4trans(l,i,k,j) = D4trans(i,j,k,l)
      D4trans(l,k,i,j) = D4trans(i,j,k,l)
      D4trans(l,j,k,i) = D4trans(i,j,k,l)
      End Do
      End Do
      End Do
      End Do
C!
C!---- Assign transformed values.
c       x = q
      call dcopy_(nvar,q,1,x,1)
c       grad = tempgrad
      call dcopy_(nvar,tempgrad,1,grad,1)

c       Hess = Temp
      call dcopy_(nvar*nvar,temp,1,Hess,1)

c       D3 = D3trans
      call dcopy_(nvar*nvar*nvar,D3trans,1,D3,1)

c       D4 = D4trans
      call dcopy_(nvar*nvar*nvar*nvar,D4trans,1,D4,1)
C!
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(var)
         Call Unused_integer(max_term)
      End If
      End
