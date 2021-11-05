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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************
!!-----------------------------------------------------------------------!
!!
      Subroutine funcval(x,coef,ipow,fval,nterm,nvar)
!!
!!  Purpose:
!!    Return function value at x.
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
      Real*8 x (nvar)
      Real*8  coef( nterm)
      Integer ipow (nterm,nvar)
      Real*8  sum,prod,fval
!!
      sum = 0.0d0
      Do iterm = 1,nterm
      prod = 1.0d0
      Do ivar = 1,nvar
      nsum = ipow(iterm,ivar)
      prod = prod*(x(ivar)**(nsum))
      End Do
      sum = sum+coef(iterm)*prod
      End Do
      fval = sum
!!
      End
!!
!!-----------------------------------------------------------------------!


!!-----------------------------------------------------------------------!
!!
      Subroutine gradient(x,coef,ipow,grad,nterm,nvar)
!!
!!  Purpose:
!!    Return gradient at x.
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
      Real*8 x (nvar)
      Real*8 coef (nvar)
      Integer ipow (nterm,nvar)
      Real*8 grad ( nvar)
      Real*8    sum,prod
      Logical   ijEq
      Real*8    rfactor
!!
      Do ivar = 1,nvar
      sum = 0.0d0
      Do iterm = 1,nterm
      prod = 1.0d0
      Do jvar = 1,nvar
      ijEq = ( ivar.eq.jvar )
      If (( ipow(iterm,jvar).ge.1 ).and.ijEq ) Then
      nder = 1
      nsum = ipow(iterm,jvar)-nder
      Else If (( ipow(iterm,jvar).ge.0 ).and.                           &
     &          ( .not.ijEq )) Then
      nder = 0
      nsum = ipow(iterm,jvar)
      Else
      nder =-1
      nsum = 0
      End If
      Call factor(nsum,nder,rfactor)
      prod = prod*rfactor*(x(jvar)**(nsum))
      End Do
      sum = sum+coef(iterm)*prod
      End Do
      grad(ivar) = sum
      End Do
!!
      End
!!
!!-----------------------------------------------------------------------!


!!-----------------------------------------------------------------------!
!!
      Subroutine Hessian(x,coef,ipow,Hess,nterm,nvar)
!!
!!  Purpose:
!!    Return Hessian at x.
!!
!!  Uses:
!!    LinAlg
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
      Real*8 x (nvar)
      Real*8 coef (nvar)
      Integer ipow (nterm,nvar)
      Real*8  Hess (nvar,nvar)
      Real*8   sum,prod
      Logical  ijEq,ikEq,jkEq
      Real*8   rfactor
!!
      Do ivar = 1,nvar
      Do jvar = ivar,nvar
      ijEq = ( ivar.eq.jvar )
      sum = 0.0d0
      Do iterm = 1,nterm
      prod = 1.0d0
      Do kvar = 1,nvar
      ikEq = ( ivar.eq.kvar )
      jkEq = ( jvar.eq.kvar )
      If (( ipow(iterm,kvar).ge.2 ).and.                                &
     &             ikEq.and.jkEq ) Then
      nder = 2
      nsum = ipow(iterm,kvar)-nder
      Else If (( ipow(iterm,kvar).ge.1 ).and.                           &
     &                        ( ikEq.or.jkEq ).and.                     &
     &     ( .not.ijEq )) Then
      nder = 1
      nsum = ipow(iterm,kvar)-nder
      Else If (( ipow(iterm,kvar).ge.0 ).and.                           &
     &                        ( .not.ikEq ).and.( .not.jkEq )) Then
      nder = 0
      nsum = ipow(iterm,kvar)
      Else
      nder =-1
      nsum = 0
      End If
      Call factor(nsum,nder,rfactor)
      prod = prod*rfactor*(x(kvar)**(nsum))
      End Do
      sum = sum+coef(iterm)*prod
      End Do
      Hess(ivar,jvar) = sum
      Hess(jvar,ivar) = sum
      End Do
      End Do
!!
      End
!!
!!-----------------------------------------------------------------------!


!!-----------------------------------------------------------------------!
!!
      Subroutine thirdDer(x,coef,ipow,D3,nterm,nvar)
!!
!!  Purpose:
!!    Return third derivatives at x.
!!
!!  Uses:
!!    LinAlg
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
      Real*8 x (nvar)
      Real*8 coef (nvar)
      Integer ipow (nterm,nvar)
      Real*8 D3 (nvar,nvar,nvar)
      Real*8   sum,prod
      Logical  ilEq,jlEq,klEq
      Real*8   rfactor
!!
!       D3 = 0.0d0
      call dcopy_(nvar*nvar*nvar,[0.0d0],0,D3,1)
      Do ivar = 1,nvar
      Do jvar = ivar,nvar
      Do kvar = jvar,nvar
      sum = 0.0d0
      Do iterm = 1,nterm
      prod = 1.0d0
      Do lvar = 1,nvar
      ilEq = ( ivar.eq.lvar )
      jlEq = ( jvar.eq.lvar )
      klEq = ( kvar.eq.lvar )
      If (( ipow(iterm,lvar).ge.3 ).and.                                &
     &                            ilEq.and.jlEq.and.klEq) Then
      nder = 3
      nsum = ipow(iterm,lvar)-nder
      Else If (( ipow(iterm,lvar).ge.2 ).and.                           &
     &                 (( ilEq.and.jlEq.and.( .not.klEq )).or.          &
     &                  ( jlEq.and.klEq.and.( .not.ilEq )))) Then
      nder = 2
      nsum = ipow(iterm,lvar)-nder
      Else If (( ipow(iterm,lvar).ge.1 ).and.                           &
     &                 (( ilEq.and.(( .not.jlEq ).and.                  &
     &                  ( .not.klEq ))).or.                             &
     &                  ( jlEq.and.(( .not.ilEq ).and.                  &
     &                  ( .not.klEq ))).or.                             &
     &                  ( klEq.and.(( .not.ilEq ).and.                  &
     &                  ( .not.jlEq ))))) Then
      nder = 1
      nsum = ipow(iterm,lvar)-nder
      Else If (( ipow(iterm,lvar).ge.0 ).and.                           &
     &                             ( .not.ilEq ).and.                   &
     &     ( .not.jlEq ).and.( .not.klEq )) Then
      nder = 0
      nsum = ipow(iterm,lvar)
      Else
      nder =-1
      nsum = 0
      End If
      Call factor(nsum,nder,rfactor)
      prod = prod*rfactor*(x(lvar)**(nsum))
      End Do
      sum = sum+coef(iterm)*prod
      End Do
      D3(ivar,jvar,kvar) = sum
      D3(kvar,ivar,jvar) = sum
      D3(jvar,kvar,ivar) = sum
      D3(jvar,ivar,kvar) = sum
      D3(kvar,jvar,ivar) = sum
      D3(ivar,kvar,jvar) = sum
      End Do
      End Do
      End Do
!!
      End
!!
!!-----------------------------------------------------------------------!


!!-----------------------------------------------------------------------!
!!
      Subroutine fourthDer(x,coef,ipow,D4,nterm,nvar)
!!
!!  Purpose:
!!    Return fourth derivatives at x.
!!
!!  Uses:
!!    LinAlg
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
      Real*8 x ( nvar)
      Real*8 coef ( nvar)
      Integer ipow (nterm,nvar)
      Real*8 D4 ( nvar, nvar, nvar, nvar)
      Real*8  sum,prod
      Logical  imEq,jmEq,kmEq,lmEq
      Real*8   rfactor
!!
!       D4 = 0.0d0
      call dcopy_(nvar*nvar*nvar*nvar,[0.0d0],0,D4,1)
      Do ivar = 1,nvar
      Do jvar = ivar,nvar
      Do kvar = jvar,nvar
      Do lvar = kvar,nvar
      sum = 0.0d0
      Do iterm = 1,nterm
      prod = 1.0d0
      Do mvar = 1,nvar
      imEq = ( ivar.eq.mvar )
      jmEq = ( jvar.eq.mvar )
      kmEq = ( kvar.eq.mvar )
      lmEq = ( lvar.eq.mvar )
      If (( ipow(iterm,mvar).ge.4 ).and.                                &
     &                               imEq.and.jmEq.and.                 &
     &                               kmEq.and.lmEq ) Then
      nder = 4
      nsum = ipow(iterm,mvar)-nder
      Else If (( ipow(iterm,mvar).ge.3 ).and.                           &
     &                              (( imEq.and.jmEq.and.               &
     &                           kmEq.and.( .not.lmEq )).or.            &
     &                               ( jmEq.and.kmEq.and.               &
     &                           lmEq.and.( .not.imEq )))) Then
      nder = 3
      nsum = ipow(iterm,mvar)-nder
      Else If (( ipow(iterm,mvar).ge.2 ).and.                           &
     &        (( imEq.and.jmEq.and.( .not.kmEq ).and.                   &
     &           ( .not.lmEq )).or.                                     &
     &         ( jmEq.and.kmEq.and.( .not.imEq ).and.                   &
     &           ( .not.lmEq )).or.                                     &
     &         ( kmEq.and.lmEq.and.( .not.imEq ).and.                   &
     &     ( .not.jmEq )))) Then
      nder = 2
      nsum = ipow(iterm,mvar)-nder
      Else If (( ipow(iterm,mvar).ge.1 ).and.                           &
     &       (( imEq.and.(( .not.jmEq ).and.( .not.kmEq ).and.          &
     &     ( .not.lmEq ))).or.                                          &
     &        ( jmEq.and.(( .not.imEq ).and.( .not.kmEq ).and.          &
     &     ( .not.lmEq ))).or.                                          &
     &        ( kmEq.and.(( .not.imEq ).and.( .not.jmEq ).and.          &
     &     ( .not.lmEq ))).or.                                          &
     &        ( lmEq.and.(( .not.imEq ).and.( .not.jmEq ).and.          &
     &     ( .not.kmEq ))))) Then
      nder = 1
      nsum = ipow(iterm,mvar)-nder
      Else If (( ipow(iterm,mvar).ge.0 ).and.                           &
     &                              ( .not.imEq ).and.                  &
     &     ( .not.jmEq ).and.( .not.kmEq ).and.( .not.lmEq )) Then
      nder = 0
      nsum = ipow(iterm,mvar)
      Else
      nder =-1
      nsum = 0
      End If
      Call factor(nsum,nder,rfactor)
      prod = prod*rfactor*(x(mvar)**(nsum))
      End Do
      sum = sum+coef(iterm)*prod
      End Do
      D4(ivar,jvar,kvar,lvar) = sum
      D4(ivar,kvar,jvar,lvar) = sum
      D4(ivar,jvar,lvar,kvar) = sum
      D4(ivar,lvar,kvar,jvar) = sum
      D4(ivar,lvar,jvar,kvar) = sum
      D4(ivar,kvar,lvar,jvar) = sum
!!
      D4(jvar,ivar,kvar,lvar) = sum
      D4(jvar,kvar,ivar,lvar) = sum
      D4(jvar,ivar,lvar,kvar) = sum
      D4(jvar,lvar,kvar,ivar) = sum
      D4(jvar,lvar,ivar,kvar) = sum
      D4(jvar,kvar,lvar,ivar) = sum
!!
      D4(kvar,ivar,jvar,lvar) = sum
      D4(kvar,jvar,ivar,lvar) = sum
      D4(kvar,ivar,lvar,jvar) = sum
      D4(kvar,lvar,jvar,ivar) = sum
      D4(kvar,lvar,ivar,jvar) = sum
      D4(kvar,jvar,lvar,ivar) = sum
!!
      D4(lvar,ivar,jvar,kvar) = sum
      D4(lvar,jvar,ivar,kvar) = sum
      D4(lvar,kvar,jvar,ivar) = sum
      D4(lvar,ivar,kvar,jvar) = sum
      D4(lvar,kvar,ivar,jvar) = sum
      D4(lvar,jvar,kvar,ivar) = sum
      End Do
      End Do
      End Do
      End Do
!!
      End
