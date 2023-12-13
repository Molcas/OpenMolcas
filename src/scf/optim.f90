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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine Optim(E_Pred,G,H,C,n,nDim)
!***********************************************************************
!                                                                      *
!    purpose: Solve a set of non-linear equations to get interpolation *
!             coefficients for damping                                 *
!                                                                      *
!     input:                                                           *
!       G       : linear terms                                         *
!       H       : bilinear terms                                       *
!       n       : size of matrices                                     *
!       nDim    : leading dimension for H                              *
!                                                                      *
!     output:                                                          *
!       C       : interpolation coefficients                           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! The KISS method is employed here. Do a scan of the energy surface    *
! along lines parallel to the edges. Reduce the steplength gradually.  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark                                                     *
!     University of Lund, Sweden, 2003                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************
      use Constants, only: Zero, Half, One
      Implicit None
!----------------------------------------------------------------------*
! Dummy arguments.                                                     *
!----------------------------------------------------------------------*
      Integer n,nDim
      Real*8  G(nDim),H(nDim,nDim),C(nDim)
!----------------------------------------------------------------------*
! Local variables.                                                     *
!----------------------------------------------------------------------*
      Real*8  Step,Eref,Eplus,Eminus,Step_m,Step_p,E_Pred
      Real*8  Optim_E, Ci, Cj
      Real*8  sum,fact
      Logical DidChange
      Integer Iter
      Integer i,j
      Logical Debug
      Logical Debug2
!----------------------------------------------------------------------*
! Initialize                                                           *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
      Debug=.true.
      Debug2=.True.
#else
      Debug=.false.
      Debug2=.false.
#endif
!----------------------------------------------------------------------*
! Make first guess.                                                    *
!----------------------------------------------------------------------*
      Do i=1,n
         C(i)=Zero
      End Do
      j=1
      Do i=1,n
         If( G(i)+Half*H(i,i) .lt. G(j)+half*H(j,j) ) j=i
      End Do
      C(j)=0.9d0
      Do i=1,n
         If(i.ne.j) C(i)=(One-C(j))/(n-1)
      End Do
      if(Debug) then
        Write(6,'(a)') 'Start C:'
        Write(6,'(6F15.6)') (C(i),i=1,n)
        Write(6,'(a)') 'Start G:'
        Write(6,'(6F26.16)') (G(i),i=1,n)
        Call RecPrt('H',' ',H,nDim,nDim)
      endif
!----------------------------------------------------------------------*
! Compute start energy.                                                *
!----------------------------------------------------------------------*
      Eref=Zero
      Do i=1,n
         Eref=Eref+C(i)*G(i)
         Do j=1,n
            Eref=Eref+Half*C(i)*C(j)*H(i,j)
         End Do
      End Do
      if(Debug) then
        Write(6,'(a,F24.16)') 'Eref = ',Eref
      endif
!----------------------------------------------------------------------*
! Do a scan in all pairwise directions                                 *
!----------------------------------------------------------------------*
      Iter=0
      DidChange=.true.
      Step=0.10d0
      Eminus=One
      Eplus =One
!
      Do While(Iter.lt.500 .and. DidChange)
!
         Iter=Iter+1
         DidChange=.false.
         Do i=1,n-1
            Do j=i+1,n
!
               Ci = C(i)
               Cj = C(j)
               Step_p = Min(Step,One-C(i),C(j))
               C(i)=C(i)+Step_p
               C(j)=C(j)-Step_p
!
               Eplus=optim_E(C,G,H,n,nDim)
!
               C(i)=Ci
               C(j)=Cj
!
               Step_m = Min(Step,C(i),One-C(j))
!
               C(i)=C(i)-Step_m
               C(j)=C(j)+Step_m
!
               EMinus=optim_E(C,G,H,n,nDim)
!
               C(i)=Ci
               C(j)=Cj
!
       if(Debug2) then
               Write(6,'(a,2I3,3F26.16)') 'i,j,Eref,Eplus,Eminus',i,j,Eref,Eplus,Eminus
               Write (6,'(a,2F26.16)') 'Step_p, Step_m=',Step_p,Step_m
               Write (6,*) 'Eplus.lt.EMinus=',Eplus.lt.EMinus
               Write (6,*) 'Eplus.lt.Eref=',Eplus.lt.Eref
               Write (6,*) 'Eminus.lt.Eref=',Eminus.lt.Eref
       endif
               If (Abs(Eplus-EMinus).gt.1.0D-12) Then
               If(Eplus.lt.Eminus) Then
                  If(Eplus.lt.Eref) Then
                     C(i)=C(i)+Step_p
                     C(j)=C(j)-Step_p
                     Eref=Eplus
                     DidChange=.true.
                  End If
               Else
                  If(Eminus.lt.Eref) Then
                     C(i)=C(i)-Step_m
                     C(j)=C(j)+Step_m
                     Eref=Eminus
                     DidChange=.true.
                  End If
               End If
               End If

            End Do ! j
         End Do ! i

         If(.not.DidChange) Then
            If(Step.gt.0.9d-4) Then
               Step=0.1d0*Step
               DidChange=.true.
               If (Debug2) Then
                  Write(6,*) 'Step is',Step
               End If
            End If
         End If
!
!        Check that the constraint is fullfilled.
!
         sum=Zero
         Do i=1,n
            If(C(i).lt.Zero) C(i)=Zero
            If(C(i).gt.One) C(i)=One
            sum=sum+C(i)
         End Do
#ifdef _DEBUGPRINT_
         Write(6,*) 'optim: sum-1',sum-One
#endif
         fact=One/sum
         Do i=1,n
            C(i)=fact*C(i)
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*) 'ERef=',ERef
#endif
      E_Pred=ERef

      End subroutine Optim


      Real*8 Function optim_E(C,G,H,n,nDim)
      Use Constants, only: Zero, Half
      Implicit None
      Integer n, nDim
      Real*8 C(nDim), G(nDim), H(nDim,nDim)

      Integer k, m
      Real*8 Tmp
!
      Optim_E=Zero
      Do k=1,n
         Tmp = Zero
         Do m=1,n
            Tmp=Tmp+Half*(C(k)*C(m)*H(k,m))
         End Do
         Optim_E=Optim_E+C(k)*G(k)+Tmp
      End Do
      End Function optim_E
