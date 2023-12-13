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
! Copyright (C) 1994, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!               2018, Ignacio Fdez. Galvan                             *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine SOrUpV(V,lvec,W,Mode,UpTp)
!***********************************************************************
!     for Ref., see T.H. Fischer and J. Almloef, JPC 96, 9768 (1992)   *
!               doi:10.1021/j100203a036                                *
!                                                                      *
!     purpose: Second Order Updated vector V using BFGS and            *
!              diagonal Hessian of Orbital Rotations                   *
!     input  : V          -> input vector                              *
!              HDiag      -> initial diagonal Hessian                  *
!              lvec       -> lengths of vectors delta, grad, HDiag & V *
!              Mode       -> update mode, see below                    *
!     The routine expects, that grad(n) and delta(n-1) are already     *
!     stored on the appropriate linked lists LLdGrd & LLDelt.          *
!     output:  W          -> H(n)V with H(n) SOrUp (inverse) Hessian   *
!                            after n-1 updates                         *
!                            mem has to be allocated from caller side  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Adapted to work both for d=-H(-1)g and g = -Hd                   *
!     depending on the value of MODE                                   *
!     Mode='DISP':  d=-H(-1)g    BFGS update for the inverse Hessian   *
!     Mode='GRAD':  g=-Hd        DFP  update for the Hessian           *
!     2017, Roland Lindh, Harvard, Cambridge                           *
!                                                                      *
!     Adapted to do both BFGS and DFP updates depending on the value   *
!     of UPTP                                                          *
!     UpTp='BFGS':  BFGS update                                        *
!     UpTp='DFP ':  DFP update                                         *
!     Mode='DISP':  d=-H(-1)g    update for the Hessian                *
!     Mode='GRAD':  g=-Hd        update for the inverse Hessian        *
!***********************************************************************
      use LnkLst, only: SCF_V
      use LnkLst, only: LLdGrd,LLDelt,LLy
!     only tentatively this Module
      use InfSO, only: IterSO
      use InfSCF, only: Iter, TimFld
      use SCF_Arrays, only: HDiag
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      use Files, only: LuDel, LuDGd
      Implicit None
!
!     declaration subroutine parameters
      Integer lvec
      Real*8 W(lVec), V(lvec)
      Character*4 Mode, UpTp
!
!     declaration local variables
      Character(LEN=8), Save:: Mode_Old
      Integer ipdel,ipdgd,ipynm1
      Integer i,it,inode,leny
      Logical updy
      Real*8 S(6),T(4)

!     declarations of functions
      Integer, External:: LstPtr, LLLen
      real*8, External:: ddot_

      Real*8 Cpu1,Cpu2,Tim1,Tim2,Tim3
      Real*8 :: Thr=1.0D-9
      Integer LL1, LL2, Lu1
      Real*8, Dimension(:), Allocatable:: SOGrd, SODel, SOScr
      Logical Inverse_H

      Interface
         SubRoutine yHx(X,Y,nXY)
         Integer nXY
         Real*8, Target:: X(nXY), Y(nXY)
         End SubRoutine yHx
      End Interface
!
      If (DDot_(lVec,V,1,V,1)==Zero) Then
         W(:)=Zero
         Return
      End If
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
      S(:)=Zero
      T(:)=Zero
!
!     Dummy initializations
!
      Inverse_H=.False.
      Lu1=-1
      LL1=-1
      LL2=-1
!
!     This section will control the mode
!
      If (Mode.eq.'DISP') Then
!
!        Mode for computation d = -H(-1)g
!
         Inverse_H=.True.
         Lu1=LuDel
         LL1=LLDelt
         LL2=LLdGrd
!
      Else If (Mode.eq.'GRAD') Then
!
!        Mode for computation g = -Hd
!
         Inverse_H=.False.
         Lu1=LudGd
         LL1=LLdGrd
         LL2=LLDelt
!
      End If
!
      If (Lu1.lt.0) Then
         Write (6,*) 'SOrUpV: Illegal mode'
         Call Abend()
      End If
!
      If ((UpTp.ne.'BFGS').and.(UpTp.ne.'DFP ')) Then
         Write (6,*) 'SOrUpV: Illegal update type',UpTp
         Call Abend()
      End If
!
      If (iterso.gt.1.AND.Mode//UpTp.ne.Mode_Old) Then
         Write (6,*) 'IterSO=',IterSO
         Write (6,*) 'Mode_Old:',Mode_Old
         Write (6,*) 'Mode//UpTp:',Mode//UpTp
         Write (6,*) 'SOrUpV: Illegal mode switch'
         Call Abend()
      End If
!
!     Steps denoted as in the reference paper!
!
!     (1): initialize w=HDiag*v
!
      If (Inverse_H) Then
         Do i=1,lvec
            If (Abs(HDiag(i)).lt.Thr) Then
               W(i)=1.0D2*V(i)
            Else
               W(i)=V(i)/HDiag(i)
            End If
         End Do
      Else
!        Do i=1,lvec
!           W(i)=HDiag(i)*V(i)
!        End Do
         Call yHx(V,W,lvec)
      End If
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*)
      Call Check_Vec(W,Size(W),'H_{n-1}v')
      Call NrmClc(V,lVec,'SORUPV','V')
      Call NrmClc(HDiag,lVec,'SORUPV','HDiag')
      Call NrmClc(W,lVec,'SORUPV','W')
#endif

!#define _DIAGONAL_ONLY_
#ifdef _DIAGONAL_ONLY_
      If (.True.) Then
         Write (6,*) ' SorUpV: Only diagonal approximation'
#else
      If (iterso.eq.1) Then
#endif
!        (2): no further work required this time
         Call Timing(Cpu2,Tim1,Tim2,Tim3)
         TimFld( 7) = TimFld( 7) + (Cpu2 - Cpu1)
         Mode_Old=Mode//UpTp
         Return
      End If
!
!     need three vectors of length lvec as scratch space
!
      Call mma_allocate(SOGrd,lvec,Label='SOGrd')
      Call mma_allocate(SODel,lvec,Label='SODel')
      Call mma_allocate(SOScr,lvec,Label='SOScr')
!
      Call GetNod(iter-1,LL2,inode)
      If (inode.eq.0) Call Error_Handling()
      Call iVPtr(SOGrd,lvec,inode)
!
!     (3b): initialize y(n-1)=HDiag*dGrd(n-1) ...
!
      If (Inverse_H) Then
        Do i=1,lvec
           If (Abs(HDiag(i)).lt.Thr) Then
              SOScr(i)=1.0D2*SOGrd(i)
           Else
              SOScr(i)=SOGrd(i)/HDiag(i)
           End If
        End Do
      Else
!       Do i=1,lvec
!          SOScr(i)=HDiag(i)*SOGrd(i)
!       End Do
        Call yHx(SOGrd,SOScr,lvec)
      End If
#ifdef _DEBUGPRINT_
!     Call RecPrt('Init y(n-1)',' ',SOScr,1,lVec)
      Call NrmClc(SOScr,lVec,'SOrUpV','Init y(n-1)')
      Call Check_vec(SOScr,lVec,'y_{n-1}=H_{0}Delta_{n-1}')
#endif
!
!     and store it on appropriate linked list
!
      leny=LLLen(LLy)
      Call PutVec(SOScr,lvec,iter-1,'NOOP',LLy)
      If (leny.eq.LLLen(LLy)) Then
!        already there, so we don't have to recalculate later
         updy=.False.
      Else
!        new listhead generated in ipy, we have to update
         updy=.True.
      End If
!
!     (4): now loop over 1..n-2 iterations.
!
#ifdef _DEBUGPRINT_
      Write (6,*) 'IterSO=',IterSO
#endif
      Do it=iter-iterso+1,iter-2
!
!        fetch delta(i), dGrd(i) and y(i) from corresponding LLists
!
         Call GetNod(it,LL1,inode)
         If (inode.eq.0) Call Error_Handling()
         Call iVPtr(SODel,lvec,inode)
!
         Call GetNod(it,LL2,inode)
         If (inode.eq.0) Call Error_Handling()
         Call iVPtr(SOGrd,lvec,inode)
!
         Call GetNod(it,LLy,inode)
         If (inode.eq.0) Call Error_Handling()
         Call iVPtr(SOScr,lvec,inode)
!
!        calculate S_k and T_k dot products.
!        (note that S(2) is the inverse of the one in the paper
!
         S(1)=ddot_(lvec,SODel,1,SOGrd,1)
         If (Abs(S(1)).lt.Thr) Then
            S(1)=Zero
            !S(1)=One/Thr
         Else
            S(1)=One/S(1)
         End If
         S(2)=ddot_(lvec,SOGrd,1,SOScr,1)
         S(3)=ddot_(lvec,SODel,1,V,1)
         S(4)=ddot_(lvec,SOScr,1,V,1)
         If (updy) Then
!
!           here we have to reload dGrd(n-1) from llist, but this
!           for sure is a memory hit, since it was put there last
!
            ipdgd=LstPtr(iter-1,LL2)
!
            S(5)=ddot_(lvec,SODel,1,SCF_V(ipdgd)%A,1)
            S(6)=ddot_(lvec,SOScr,1,SCF_V(ipdgd)%A,1)
         Else
            S(5)=Zero
            S(6)=Zero
         End If
#ifdef _DEBUGPRINT_
         Write (6,*) 'it=',it
         Write (6,*) '(S(i),i=1,6)=',(S(i),i=1,6)
#endif
!
         If ((Mode_Old.eq.'DISPBFGS').or.(Mode_Old.eq.'GRADDFP ')) Then
            If (Abs(S(2)).lt.Thr) Then
               S(2)=One
               !S(2)=One+S(1)/Thr
            Else
               S(2)=One+S(1)*S(2)
            End If
            T(2)=S(1)*S(3)
            T(4)=S(1)*S(5)
            T(1)=S(2)*T(2)-S(1)*S(4)
            T(3)=S(2)*T(4)-S(1)*S(6)
         Else If ((Mode_Old.eq.'DISPDFP ').or.(Mode_Old.eq.'GRADBFGS')) Then
!           DFP update (or BFGS for the Hessian) use
!           a slightly different and simpler formula
            If (Abs(S(2)).lt.Thr) Then
               S(2)=Zero
               S(2)=One/Thr
            Else
               S(2)=One/S(2)
            End If
            T(1)=S(1)*S(3)
            T(2)=S(2)*S(4)
            T(3)=S(1)*S(5)
            T(4)=S(2)*S(6)
         End If
!
!        Compute w and y(n-1)
!
#ifdef _DEBUGPRINT_
         Write (6,*) '(T(i),i=1,4)=',(T(i),i=1,4)
         Call Check_Vec(W,Size(W),'W(0)')
         Call daxpy_(lvec, T(1),SODel,1,W,1)
         Call Check_Vec(W,Size(W),'W(1)')
         Call daxpy_(lvec,-T(2),SOScr,1,W,1)
         Call Check_Vec(W,Size(W),'W(2)')
#else
         W(:)=W(:)+T(1)*SODel(:)-T(2)*SOScr(:)
#endif
         If (updy) Then
!
!           here we have to reload y(n-1) from llist, but this
!           for sure is a memory hit, since it was put there last
!           -> we operate directly on the memory cells of the LList,
!           where y(n-1) resides.
!
            ipynm1=LstPtr(iter-1,LLy)
!           Call daxpy_(lvec, T(3),SODel,1,SCF_V(ipynm1)%A,1)
!           Call daxpy_(lvec,-T(4),SOScr,1,SCF_V(ipynm1)%A,1)
            SCF_V(ipynm1)%A(:) = SCF_V(ipynm1)%A(:) + T(3)*SODel(:) - T(4)*SOScr(:)
         End If
!
      End Do
!
!     (5): reload y(n-1), delta(n-1) & dGrd(n-1) from linked list.
!     these all are memory hits, of course
!
      ipynm1=LstPtr(iter-1,LLy)
      ipdel =LstPtr(iter-1,LL1)
      ipdgd =LstPtr(iter-1,LL2)
#ifdef _DEBUGPRINT_
      Write (6,*)
!     Call RecPrt('y(n-1)',' ',SCF_V(ipynm1)%A,1,lVec)
      Call NrmClc(SCF_V(ipynm1)%A,lVec,'SOrUpV','y(n-1)')
      If (Mode.eq.'DISP') Then
!        Call RecPrt('dX(n-1)',' ',SCF_V(ipdel)%A,1,lVec)
!        Call RecPrt('dg(n-1)',' ',SCF_V(ipdgd)%A,1,lVec)
         Call NrmClc(SCF_V(ipdel)%A,lVec,'SOrUpV','dX(n-1)')
         Call NrmClc(SCF_V(ipdgd)%A,lVec,'SOrUpV','dg(n-1)')
      Else
!        Call RecPrt('dg(n-1)',' ',SCF_V(ipdel)%A,1,lVec)
!        Call RecPrt('dX(n-1)',' ',SCF_V(ipdgd)%A,1,lVec)
         Call NrmClc(SCF_V(ipdel)%A,lVec,'SOrUpV','dg(n-1)')
         Call NrmClc(SCF_V(ipdgd)%A,lVec,'SOrUpV','dX(n-1)')
      End If
      Write (6,*)
      Call Check_Vec(SCF_V(ipdel)%A,lvec,'delta_{n-1}')
      Call Check_Vec(SCF_V(ipdgd)%A,lvec,'Delta_{n-1}')
#endif
!
!     calculate diverse dot products...
!
      S(1)=ddot_(lvec,SCF_V(ipdel)%A,1,SCF_V(ipdgd)%A,1)
!     Write (6,*) 'S(1)=',S(1)
      If (Abs(S(1)).lt.Thr) Then
         S(1)=Zero
         !S(1)=One/Thr
      Else
         S(1)=One/S(1)
      End If
!     Write (6,*) 'S(1)=Alpha=',S(1)
!     Call Check_vec(SCF_V(ipynm1)%A,lVec,'y_{n-1}=H_{n-1}Delta_{n-1}')
      S(2)=ddot_(lvec,SCF_V(ipdgd)%A,1,SCF_V(ipynm1)%A,1)
!     Write (6,*) 'S(2)=Delta_{n-1}^T*y_{n-1}=',S(2)
      S(3)=ddot_(lvec,SCF_V(ipdel)%A,1,V,1)
!     Write (6,*) 'S(3)=delta_{n-1}^T*v=',S(3)
      S(4)=ddot_(lvec,SCF_V(ipynm1)%A,1,V,1)
!     Write (6,*) 'S(4)=y_{n-1}^T*v=',S(4)

#ifdef _DEBUGPRINT_
      Write (6,*) '(S(i),i=1,4)=',(S(i),i=1,4)
#endif

      If ((Mode_Old.eq.'DISPBFGS').or.(Mode_Old.eq.'GRADDFP ')) Then
         If (Abs(S(2)).lt.Thr) Then
            S(2)=One
            !S(2)=One+S(1)/Thr
         Else
            S(2)=One+S(1)*S(2)
         End If
!        Write (6,*) 'S(2)(1+...)=',S(2)
!        Write (6,*)
         T(2)=S(1)*S(3)
         T(1)=S(2)*T(2)-S(1)*S(4)
!        Write (6,*) 'T(:)=',T(:)
      Else If ((Mode_Old.eq.'DISPDFP ').or.(Mode_Old.eq.'GRADBFGS')) Then
         If (Abs(S(2)).lt.Thr) Then
            S(2)=Zero
            S(2)=One/Thr
         Else
            S(2)=One/S(2)
         End If
         T(1)=S(1)*S(3)
         T(2)=S(2)*S(4)
      End If
!
!     update the vector w
!
#ifdef _DEBUGPRINT_
      Write (6,*) '(T(i),i=1,2)=',(T(i),i=1,2)
      Call Check_Vec(W,Size(W),'W(2), again')
      Call daxpy_(lvec, T(1),SCF_V(ipdel)%A,1,W,1)
      Call Check_Vec(SCF_V(ipdel)%A,lvec,'delta_{n-1}')
      Call Check_Vec(W,Size(W),'W(3)')
      Call daxpy_(lvec,-T(2),SCF_V(ipynm1)%A,1,W,1)
      Call Check_Vec(SCF_V(ipynm1)%A,lvec,'y_{n-1}')
!     Call RecPrt('The final W array',' ',W,1,lVec)
      Call NrmClc(W,lVec,'SOrUpV','The final W array')
      Call Check_Vec(W,Size(W),'W(final)')
#else
      W(:)=W(:)+T(1)*SCF_V(ipdel)%A(:)-T(2)*SCF_V(ipynm1)%A(:)
#endif
!
!     so, we've made it, let's clean up workbench
!
      Call mma_deallocate(SOGrd)
      Call mma_deallocate(SODel)
      Call mma_deallocate(SOScr)
!
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 7) = TimFld( 7) + (Cpu2 - Cpu1)
      Return
!
!-----Error handling
!
!     Hmmm, no entry found in LList, that's strange
      Contains
      Subroutine Error_handling()
      Write (6,*) 'SOrUpV: no entry found in LList'
      Call Abend()
      End Subroutine Error_handling
#ifdef _DEBUGPRINT_
      Subroutine Check_Vec(Vec,nVec,Label)
      Integer nVec
      Real*8 Vec(nVec), DD
      Real*8, External:: DDot_
      Character(Len=*) Label
      DD=Sqrt(DDot_(nVec,Vec,1,Vec,1))
      Write (6,*) 'Norm of ',Label,' is : ',DD
      End Subroutine Check_Vec
#endif

      End Subroutine Sorupv
