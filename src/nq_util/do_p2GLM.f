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
      Subroutine Do_P2GLM(P2mo,np2act,D1mo,nd1mo,TabMO,mAO,mGrid,nMOs,
     &                   P2_ontop,nP2_ontop,RhoI,RhoA,mRho,Do_Grad)
************************************************************************
*                                                                      *
* Object: Calculation P2 ontop density and its derivatives             *
*                                                                      *
* Called from: Do_batch                                                *
*                                                                      *
* Calling    : FZero                                                   *
*                                                                      *
*    INPUT:                                                            *
*   P2mo     = two-body density matrix in MO basis                     *
*   np2act   = size of P2mo                                            *
*   D1mo     = one-body density matrix in MO basis                     *
*   nd1mo    = size of D1mo                                            *
*   TabMO    = MO values computed on grid                              *
*   nMOs     = number of MO basis                                      *
*   mAO      = number of derivatives of AO...                          *
*   mGrid    = number of grid points                                   *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "nq_info.fh"
#include "real.fh"
      Real*8 P2mo(np2act),D1mo(nd1mo),TabMO(mAO,mGrid,nMOs),
     &       P2_ontop(nP2_ontop,mGrid)
      Integer iOff_Ash(0:7), iOff_Bas(0:7), iOff_BasAct(0:7)
      Real*8 RhoI(mRho,mGrid)
      Real*8 RhoA(mRho,mGrid)
      Logical Do_Grad
************************************************************************
*                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
c          write(6,*) 'nP2_ontop =', nP2_ontop
c          write(6,*) 'np2act =', np2act
c          write(6,*) 'mAO       =', mAO
c          write(6,*) 'mRho      =', mRho
      If (nP2_ontop.eq.4) Then
         If (mAO.ne.4.or.mRho.ne.4) Then
           Call WarningMessage(2,' Somthings wrong in dim. in p2cs')
           Call Abend()
       End If
      Else If (nP2_ontop.eq.6) Then
         If (mAO.ne.10.or.mRho.ne.6) Then
           Call WarningMessage(2,' Somthings wrong in dim. in p2cs')
           Call Abend()
      End If
      End If
*
      Call FZero(P2_ontop,mGrid*nP2_ontop)
      jOffA_ = 0
      jOffB_ = 0
      Do iIrrep = 0, mIrrep-1
         iOff_Ash(iIrrep)=jOffA_
         iOff_Bas(iIrrep)=jOffB_
         iOff_BasAct(iIrrep)=jOffB_ + nIsh(iIrrep) + nFro(iIrrep)
         jOffA_=jOffA_+nAsh(iIrrep)
         jOffB_=jOffB_+mBas(iIrrep)
      End Do
************************************************************************
*                                                                      *
*   P(1,...) - P_2                                                     *
*   P(2,...), P(3,...), P(4,...) - grad P_2                            *
*   Not implemented:                                                   *
*   P(5,...) - grad^2 P_2                                              *
*   P(6,...) - additional part grad^2 P_2 for CS functional            *
*   P(5) and P(6) removed                                              *
************************************************************************

************************************************************************
*          Inactive part:                                              *
************************************************************************
      NumIsh = 0
      NumAsh = 0
      Do iIrrep=0, mIrrep-1
         NumIsh = NumIsh + nISh(iIrrep)
         NumAsh = NumAsh + nAsh(iIrrep)
      End Do
*
      Do iGrid = 1, mGrid
      Do iIrrep=0, mIrrep-1
*      Write(6,*) " Symm:",iIrrep
       Do i_=1,nISh(iIrrep) + nFro(iIrrep)
        i = iOff_Bas(iIrrep) + i_
c
c         Write(6,*)  " do_p2: Inact-Inact:", iIrrep,i,
c     &               TabMO(1,iGrid,i)
c
          RhoI(1,iGrid) = RhoI(1,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(1,iGrid,i)
c      write(6,'(A15,2I3,2G15.8)')'iGrid,i,MO,RhoI',
c     &   iGrid,i,TabMO(1,iGrid,i), RhoI(1,iGrid)
!        if (Functional_type.eq.GGA_type.or.Do_Grad) then
        if (Functional_type.eq.GGA_type) then
          RhoI(2,iGrid) = RhoI(2,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(2,iGrid,i)
          RhoI(3,iGrid) = RhoI(3,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(3,iGrid,i)
          RhoI(4,iGrid) = RhoI(4,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(4,iGrid,i)
        end if
*
       End Do         ! i_
      End Do         ! iIrrep
      End Do         ! iGrid
*
      If (NumIsh.ne.0) Then
      Do iGrid = 1, mGrid
         P2_ontop(1,iGrid) = RhoI(1,iGrid)*RhoI(1,iGrid)
c         Write(6,'(A15,I3,1G28.20)')'iGrid,P2(1)=',
c     *          iGrid,P2_ontop(1,iGrid)
*
!        if (Functional_type.eq.GGA_type.or.Do_Grad) then
        if (Functional_type.eq.GGA_type) then
            P2_ontop(2,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(2,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(2)   =',P2_ontop(2,iGrid)
            P2_ontop(3,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(3,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(3)   =',P2_ontop(3,iGrid)
            P2_ontop(4,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(4,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(4)   =',P2_ontop(4,iGrid)
         end if
        if (Functional_type.eq.LDA_type.and.Do_Grad) then
!Here I must
!1. transform the 2-body density matrix to AO

!2. Loop over effective gradients
!3. Calculate P2_ontop_d(eff_Grad,iGrid)


        end if
      End Do
      End If

************************************************************************
*          Active-Inactive part:                                       *
************************************************************************
      If (NumIsh.ne.0.and.NumAsh.ne.0) Then
       Do kIrrep = 0, mIrrep-1
         Do k_ = 1, nASh(kIrrep)
            k= k_ + iOff_BasAct(kIrrep)
            Do lIrrep = 0, mIrrep-1
               Do l_ = 1, nAsh(lIrrep)
                  l= l_ + iOff_BasAct(lIrrep)
                  kl=iTri(k_ + iOff_Ash(kIrrep) ,
     *                    l_ + iOff_Ash(lIrrep) )
                  Do iGrid = 1, mGrid
                     RhoA(1,iGrid) = RhoA(1,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(1,iGrid,l)
c       Write(6,'(A35,3I3,3G15.8)') 'iGrid,k,l,D1mo(kl),Tab(k),Tab(l)=',
c     &    iGrid,k,l,D1mo(kl),TabMO(1,iGrid,k),TabMO(1,iGrid,l)
!                     if (Functional_type.eq.GGA_type.or.Do_Grad) Then
                     if (Functional_type.eq.GGA_type) Then
                     RhoA(2,iGrid) = RhoA(2,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(2,iGrid,l)
                     RhoA(3,iGrid) = RhoA(3,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(3,iGrid,l)
C        write(6,*) 'RhoA(4,iGrid) bf =', RhoA(4,iGrid)
                     RhoA(4,iGrid) = RhoA(4,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
*       Write(6,*) 'D1mo(kl),Tab(1,k),Tab(1,l)=',
*     &    D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
                     end if
                  End Do     ! iGrid
               End Do      ! l_
            End Do       ! lIrrep
         End Do        ! k_
       End Do        ! kIrrep
*
       Do iGrid = 1, mGrid
               P2_ontop(1,iGrid) = P2_ontop(1,iGrid) +
     *                           RhoI(1,iGrid)*RhoA(1,iGrid)
c       write(6,'(A20,I3,2G15.8)')'iGrid,RhoI,RhoA=',
c     *                           iGrid,RhoI(1,iGrid),RhoA(1,iGrid)
c       Write(6,'(A15,I3,1f28.20)') 'iGrid,P2(1) =',
c     *    iGrid,P2_ontop(1,iGrid)
!        if (Functional_type.eq.GGA_type.or.Do_Grad) Then
        if (Functional_type.eq.GGA_type) Then
               P2_ontop(2,iGrid) = P2_ontop(2,iGrid) +
     *                     2.0d0*RhoI(2,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(2,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(2)   =',P2_ontop(2,iGrid)
               P2_ontop(3,iGrid) = P2_ontop(3,iGrid) +
     *                     2.0d0*RhoI(3,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(3,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(3)   =',P2_ontop(3,iGrid)
               P2_ontop(4,iGrid) = P2_ontop(4,iGrid) +
     *                     2.0d0*RhoI(4,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(4,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(4)   =',P2_ontop(4,iGrid)
        end if
       End Do ! loop over grid points
      End If ! if Inactive
************************************************************************
*
*          Active-Active part:
*
************************************************************************
      If (NumAsh.ne.0) Then
      Do iIrrep = 0, mIrrep-1
         Do jIrrep = 0, mIrrep-1
         ijIrrep=iEor(iIrrep,jIrrep)
            Do kIrrep = 0, mIrrep-1
            ijkIrrep=iEor(ijIrrep,kIrrep)
*      Write(6,*) " do_p2:Irrep" ,iIrrep,jIrrep,kIrrep,ijkIrrep
               Do k_ = 1, nASh(kIrrep)
                  k = iOff_BasAct(kIrrep)+k_
                  Do l_ = 1, nASh(ijkIrrep)
                     l=iOff_BasAct(ijkIrrep)+l_
                     Do i_ = 1, nASh(iIrrep)
                        i=iOff_BasAct(iIrrep)+i_
                        Do j_ = 1, nASh(jIrrep)
                           j=iOff_BasAct(jIrrep)+j_
                           kl=iTri(k_ + iOff_Ash(  kIrrep) ,
     &                             l_ + iOff_Ash(ijkIrrep))
                           ij=iTri(i_ + iOff_Ash(iIrrep) ,
     &                             j_ + iOff_Ash(  jIrrep))
                           ijkl  = iTri(ij,kl)
c incorrect piece of code ....
c                             Fact=0.5d0
c                             if(iIrrep.eq.jIrrep) Then
c                                If(k.eq.l) Fact=1.0d0
c                             End If
c correct piece of code ....
                             Fact=0.5d0
                             if((ij.ge.kl).and.(k.eq.l)) Fact=1.0d0
                             if((kl.ge.ij).and.(i.eq.j)) Fact=1.0d0
c*****************************************************************************
c        write(6,'(4I2,F4.1,F8.5)')  i,j,k,l,Fact,P2mo(ijkl)
                           Do iGrid = 1, mGrid
c       write(6,'(I3,F5.2,4I2,5E15.8)') iGrid,Fact,i,j,k,l,
c     &  P2mo(ijkl),
c     &   TabMO(1,igrid,i),TabMO(1,igrid,j),
c     &  TabMO(1,igrid,k),TabMO(1,igrid,l)

                             P2_ontop(1,iGrid) = P2_ontop(1,iGrid)  +
     *                                    Fact*P2mo(ijkl)*
     *                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     *                             TabMO(1,igrid,i)*TabMO(1,igrid,j)
c         Write(6,'(A,1G28.20)') 'P2(1) =',P2_ontop(1,iGrid)
                             if (Functional_type.eq.GGA_type
     &                                     ) Then
!     &                                     .or.Do_Grad) Then
                              P2_ontop(2,iGrid) = P2_ontop(2,iGrid)  +
     &                                     Fact*P2mo(ijkl)*(
     &                             TabMO(2,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(2,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(2,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(2,igrid,j)
     &                             )
C            Write(6,'(A,1f28.20)') 'P2(2)   =',P2_ontop(2,iGrid)
                              P2_ontop(3,iGrid) = P2_ontop(3,iGrid) +
     &                                    Fact*P2mo(ijkl)*(
     &                             TabMO(3,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(3,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(3,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(3,igrid,j)
     &                             )
C            Write(6,'(A,1f28.20)') 'P2(3)   =',P2_ontop(3,iGrid)
                              P2_ontop(4,iGrid) = P2_ontop(4,iGrid) +
     &                                    Fact*P2mo(ijkl)*(
     &                             TabMO(4,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(4,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(4,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(4,igrid,j)
     &                             )
C            Write(6,'(A,1f28.20)') 'P2(4)   =',P2_ontop(4,iGrid)
                              end if
*
                           End Do ! iGrid
                        End Do ! j_
                     End Do ! i_
                  End Do ! l_
               End Do ! k_
            End Do ! kIrrep
         End Do ! jIrrep
      End Do ! iIrrep
      End If
*
      RETURN
      END



      Subroutine Do_P2GLM_grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &          P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
     &          list_s,nlist_s,list_bas,Index,nIndex,
     &          P2mo,np2act,D1mo,nd1mo,TabMO,list_g,P2_ontop_d,
     &          RhoI,RhoA,mRho,nMOs,CMO,nAOs)
************************************************************************
*                                                                      *
* Object: Calculation P2 ontop density and its derivatives             *
*                                                                      *
* Called from: Do_batch                                                *
*                                                                      *
* Calling    : FZero                                                   *
*                                                                      *
*    INPUT:                                                            *
*   P2mo     = two-body density matrix in MO basis                     *
*   np2act   = size of P2mo                                            *
*   D1mo     = one-body density matrix in MO basis                     *
*   nd1mo    = size of D1mo                                            *
*   TabMO    = MO values computed on grid                              *
*   nMOs     = number of MO basis                                      *
*   mAO      = number of derivatives of AO...                          *
*   mGrid    = number of grid points                                   *
*                                                                      *
************************************************************************
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
#include "nq_info.fh"
#include "real.fh"
!Error could be TabAO...
      Real*8 P2mo(np2act),D1mo(nd1mo),TabMO(mAO,mGrid,nMOs),
     &     P2_ontop(nP2_ontop,mGrid),TabAO(nTabAO),
     &     P2_ontop_d(nGrad_Eff,mGrid),CMO(nMOs,nAOs)
!      Real*8 P2AO(np2AO)
      Integer iOff_Ash(0:7), iOff_Bas(0:7), iOff_BasAct(0:7),
     &        list_s(2,nlist_s),list_bas(2,nlist_s),Index(nIndex),
     &        list_g(3,nlist_s),ipTabAO(nlist_s),
     &        mAO,nAOs,mGrid,nP2_ontop,nGrad_Eff,np2act,nd1mo,nTabAO
      Real*8, allocatable, dimension(:,:,:) :: AO_vals
!      Integer np2AO
      Real*8 RhoI(mRho,mGrid)
      Real*8 RhoA(mRho,mGrid)
      Real*8 gf1(1:3)
      Logical Do_Grad
      integer iCB,g_eff



************************************************************************
*                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
!      write(*,*) 'CMO print'
!      do i=1,nMOs
!       write(*,*) 'MO number: ',i
!       do j=1,nAOs
!       write(*,*) j,CMO(j,i)
!       end do
!      end do


      if (.false.) then
      Do ilist_s=1,nlist_s
        write(6,*) 'INFO FOR SHELL',ilist_s
        iSkal=list_s(1,ilist_s)
        write(6,*) 'iskal',iskal
        iCmp  = iSD( 2,iSkal)
        write(6,*) 'icmp',icmp
        iBas  = iSD( 3,iSkal)
        write(6,*) 'ibas',ibas
        iBas_Eff = List_Bas(1,ilist_s)
        write(6,*) 'iBas_eff',ibas_eff
        mdci  = iSD(10,iSkal)
        iShell= iSD(11,iSkal)
        write(6,*) 'ishell',ishell
        index_i=list_bas(2,ilist_s)
        write(6,*) 'index_i',index_i
        nFunc_i=iBas*iCmp
        IndMap = Index(index_i)
      end do
      end if


      If (nP2_ontop.eq.4) Then
         If (mAO.ne.4.or.mRho.ne.4) Then
           Call WarningMessage(2,' Somthings wrong in dim. in p2cs')
           Call Abend()
       End If
      Else If (nP2_ontop.eq.6) Then
         If (mAO.ne.10.or.mRho.ne.6) Then
           Call WarningMessage(2,' Somthings wrong in dim. in p2cs')
           Call Abend()
      End If
      End If
*
!       write(*,*) 'Index: ',Index
!       write(*,*) 'ilistbas', List_bas
      Call FZero(P2_ontop,mGrid*nP2_ontop)
      jOffA_ = 0
      jOffB_ = 0
      Do iIrrep = 0, mIrrep-1
         iOff_Ash(iIrrep)=jOffA_
         iOff_Bas(iIrrep)=jOffB_
         iOff_BasAct(iIrrep)=jOffB_ + nIsh(iIrrep) + nFro(iIrrep)
         jOffA_=jOffA_+nAsh(iIrrep)
         jOffB_=jOffB_+mBas(iIrrep)
      End Do
************************************************************************
*                                                                      *
*   P(1,...) - P_2                                                     *
*   P(2,...), P(3,...), P(4,...) - grad P_2                            *
*   Not implemented:                                                   *
*   P(5,...) - grad^2 P_2                                              *
*   P(6,...) - additional part grad^2 P_2 for CS functional            *
*   P(5) and P(6) removed                                              *
*                                                                      *
************************************************************************

************************************************************************
*          Inactive part:                                              *
************************************************************************
      NumIsh = 0
      NumAsh = 0
      Do iIrrep=0, mIrrep-1
         NumIsh = NumIsh + nISh(iIrrep)
         NumAsh = NumAsh + nAsh(iIrrep)
      End Do
*
      Do iGrid = 1, mGrid
      Do iIrrep=0, mIrrep-1
*      Write(6,*) " Symm:",iIrrep
       Do i_=1,nISh(iIrrep) + nFro(iIrrep)
        i = iOff_Bas(iIrrep) + i_



          RhoI(1,iGrid) = RhoI(1,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(1,iGrid,i)
c      write(6,'(A15,2I3,2G15.8)')'iGrid,i,MO,RhoI',
c     &   iGrid,i,TabMO(1,iGrid,i), RhoI(1,iGrid)
!        if (Functional_type.eq.GGA_type.or.Do_Grad) then
!        if (Functional_type.eq.GGA_type) then
!          RhoI(2,iGrid) = RhoI(2,iGrid) +
!     *                    TabMO(1,iGrid,i)*TabMO(2,iGrid,i)
!          RhoI(3,iGrid) = RhoI(3,iGrid) +
!     *                    TabMO(1,iGrid,i)*TabMO(3,iGrid,i)
!          RhoI(4,iGrid) = RhoI(4,iGrid) +
!     *                    TabMO(1,iGrid,i)*TabMO(4,iGrid,i)
*
       End Do         ! i_
      End Do         ! iIrrep
      End Do         ! iGrid
*
      If (NumIsh.ne.0) Then
      Do iGrid = 1, mGrid
         P2_ontop(1,iGrid) = RhoI(1,iGrid)*RhoI(1,iGrid)
c         Write(6,'(A15,I3,1G28.20)')'iGrid,P2(1)=',
c     *          iGrid,P2_ontop(1,iGrid)
*
!        if (Functional_type.eq.GGA_type.or.Do_Grad) then
        if (Functional_type.eq.GGA_type) then
            P2_ontop(2,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(2,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(2)   =',P2_ontop(2,iGrid)
            P2_ontop(3,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(3,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(3)   =',P2_ontop(3,iGrid)
            P2_ontop(4,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(4,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(4)   =',P2_ontop(4,iGrid)
         end if
        if (Functional_type.eq.LDA_type.and.Do_Grad) then
!Here I must
!1. transform the 2-body density matrix to AO

!2. Loop over effective gradients
!3. Calculate P2_ontop_d(eff_Grad,iGrid)


        end if
      End Do
      End If

************************************************************************
*          Active-Inactive part:                                       *
************************************************************************
      If (NumIsh.ne.0.and.NumAsh.ne.0) Then
       Do kIrrep = 0, mIrrep-1
         Do k_ = 1, nASh(kIrrep)
            k= k_ + iOff_BasAct(kIrrep)
            Do lIrrep = 0, mIrrep-1
               Do l_ = 1, nAsh(lIrrep)
                  l= l_ + iOff_BasAct(lIrrep)
                  kl=iTri(k_ + iOff_Ash(kIrrep) ,
     *                    l_ + iOff_Ash(lIrrep) )
                  Do iGrid = 1, mGrid
                     RhoA(1,iGrid) = RhoA(1,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(1,iGrid,l)
c       Write(6,'(A35,3I3,3G15.8)') 'iGrid,k,l,D1mo(kl),Tab(k),Tab(l)=',
c     &    iGrid,k,l,D1mo(kl),TabMO(1,iGrid,k),TabMO(1,iGrid,l)
!                     if (Functional_type.eq.GGA_type.or.Do_Grad) Then
                     if (Functional_type.eq.GGA_type) Then
                     RhoA(2,iGrid) = RhoA(2,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(2,iGrid,l)
                     RhoA(3,iGrid) = RhoA(3,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(3,iGrid,l)
C        write(6,*) 'RhoA(4,iGrid) bf =', RhoA(4,iGrid)
                     RhoA(4,iGrid) = RhoA(4,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
*       Write(6,*) 'D1mo(kl),Tab(1,k),Tab(1,l)=',
*     &    D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
                     end if
                  End Do     ! iGrid
               End Do      ! l_
            End Do       ! lIrrep
         End Do        ! k_
       End Do        ! kIrrep
*
       Do iGrid = 1, mGrid
               P2_ontop(1,iGrid) = P2_ontop(1,iGrid) +
     *                           RhoI(1,iGrid)*RhoA(1,iGrid)
c       write(6,'(A20,I3,2G15.8)')'iGrid,RhoI,RhoA=',
c     *                           iGrid,RhoI(1,iGrid),RhoA(1,iGrid)
c       Write(6,'(A15,I3,1f28.20)') 'iGrid,P2(1) =',
c     *    iGrid,P2_ontop(1,iGrid)
!        if (Functional_type.eq.GGA_type.or.Do_Grad) Then
        if (Functional_type.eq.GGA_type) Then
               P2_ontop(2,iGrid) = P2_ontop(2,iGrid) +
     *                     2.0d0*RhoI(2,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(2,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(2)   =',P2_ontop(2,iGrid)
               P2_ontop(3,iGrid) = P2_ontop(3,iGrid) +
     *                     2.0d0*RhoI(3,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(3,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(3)   =',P2_ontop(3,iGrid)
               P2_ontop(4,iGrid) = P2_ontop(4,iGrid) +
     *                     2.0d0*RhoI(4,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(4,iGrid)
C            Write(6,'(A,1f28.20)') 'P2(4)   =',P2_ontop(4,iGrid)
        end if
       End Do ! loop over grid points
      End If ! if Inactive
************************************************************************
*
*          Active-Active part:
*
************************************************************************
      If (NumAsh.ne.0) Then
      Do iIrrep = 0, mIrrep-1
         Do jIrrep = 0, mIrrep-1
         ijIrrep=iEor(iIrrep,jIrrep)
            Do kIrrep = 0, mIrrep-1
            ijkIrrep=iEor(ijIrrep,kIrrep)
*      Write(6,*) " do_p2:Irrep" ,iIrrep,jIrrep,kIrrep,ijkIrrep
               Do k_ = 1, nASh(kIrrep)
                  k = iOff_BasAct(kIrrep)+k_
                  Do l_ = 1, nASh(ijkIrrep)
                     l=iOff_BasAct(ijkIrrep)+l_
                     Do i_ = 1, nASh(iIrrep)
                        i=iOff_BasAct(iIrrep)+i_
                        Do j_ = 1, nASh(jIrrep)
                           j=iOff_BasAct(jIrrep)+j_
                           kl=iTri(k_ + iOff_Ash(  kIrrep) ,
     &                             l_ + iOff_Ash(ijkIrrep))
                           ij=iTri(i_ + iOff_Ash(iIrrep) ,
     &                             j_ + iOff_Ash(  jIrrep))
                           ijkl  = iTri(ij,kl)
                           Fact=0.5d0
                           if(iIrrep.eq.jIrrep) Then
                              If(k.eq.l) Fact=1.0d0
                           End If
c        write(6,'(4I2,F4.1,F8.5)')  i,j,k,l,Fact,P2mo(ijkl)
                           Do iGrid = 1, mGrid
c       write(6,'(I3,F5.2,4I2,5E15.8)') iGrid,Fact,i,j,k,l,
c     &  P2mo(ijkl),
c     &   TabMO(1,igrid,i),TabMO(1,igrid,j),
c     &  TabMO(1,igrid,k),TabMO(1,igrid,l)

!Transform i to AO:
      gf1(1:3) = 0
      Do ilist_s=1,nlist_s
        iSkal=list_s(1,ilist_s)
        iCmp  = iSD( 2,iSkal)
        iBas  = iSD( 3,iSkal)
        iBas_Eff = List_Bas(1,ilist_s)
        mdci  = iSD(10,iSkal)
        iShell= iSD(11,iSkal)
        index_i=list_bas(2,ilist_s)
        nFunc_i=iBas*iCmp

!I don't really understand the structure of the TabAO, so I move the
!info to a new allocated array.

!         write(*,*) 'ibas_eff,icmp',iBas_Eff,iCmp

         Allocate(AO_vals(1:mAO,1:mGrid,1:iBas_Eff*iCmp))
         Call Get_AO_info(TabAO(ipTabAO(ilist_s)),mAO,mGrid,iBas_Eff,
     &         iCmp,AO_vals)
!        allocate(IndMap(1:iBas*iCmp))
!        IndMap = Index(index_i)
!        write(*,*) 'IndMap',IndMap
!      write(*,*) "ao index, mo index, value"
      Do iCB_Eff = 1, iBas_Eff*iCmp
!I think this is an acceptable way to get the AO index...
         iCB = iCB_Eff + index_i - 1
!         write(*,*) 'icb',iCB
!         if(iCB.eq.2) write(*,*) 'YES'
!         write(*,*) i,iCB,CMO(iCB,i)
         Do nGrad=1,3
           gf1(nGrad) = gf1(nGrad) +
     &                  CMO(iCB,i)*AO_vals(nGrad+1,iGrid,iCB_Eff)
!           write(*,*) iCB,i,CMO(iCB,i)
         end do
      end do
      Do nGrad=1,3
        g_eff = list_g(nGrad,ilist_s)
        P2_ontop_d(g_eff,iGrid) = P2_ontop_d(g_eff,iGrid) +
     &    Fact*P2mo(ijkl)*gf1(nGrad)*TabMO(1,igrid,j)*
     &    TabMO(1,igrid,k)*TabMO(1,iGrid,l)
      end do
         Deallocate(AO_vals)
      end do

!Transform j to AO:
      gf1(1:3) = 0
      Do ilist_s=1,nlist_s
        iSkal=list_s(1,ilist_s)
        iCmp  = iSD( 2,iSkal)
        iBas  = iSD( 3,iSkal)
        iBas_Eff = List_Bas(1,ilist_s)
        mdci  = iSD(10,iSkal)
        iShell= iSD(11,iSkal)
        index_i=list_bas(2,ilist_s)
        nFunc_i=iBas*iCmp

         Allocate(AO_vals(1:mAO,1:mGrid,1:iBas_Eff*iCmp))
         Call Get_AO_info(TabAO(ipTabAO(ilist_s)),mAO,mGrid,iBas_Eff,
     &         iCmp,AO_vals)

      Do iCB_Eff = 1, iBas_Eff*iCmp
         iCB = iCB_Eff + index_i - 1
!         iCB = Index(iCB_Eff)
         Do nGrad=1,3
           gf1(nGrad) = gf1(nGrad) +
     &                  CMO(iCB,i)*AO_vals(nGrad+1,iGrid,iCB_Eff)
         end do
      end do
      Do nGrad=1,3
        g_eff = list_g(nGrad,ilist_s)
        P2_ontop_d(g_eff,iGrid) = P2_ontop_d(g_eff,iGrid) +
     &    Fact*P2mo(ijkl)*gf1(nGrad)*TabMO(1,igrid,i)*
     &    TabMO(1,igrid,k)*TabMO(1,iGrid,l)
      end do
         Deallocate(AO_vals)
      end do

!Transform k to AO:
      gf1(1:3) = 0
      Do ilist_s=1,nlist_s
        iSkal=list_s(1,ilist_s)
        iCmp  = iSD( 2,iSkal)
        iBas  = iSD( 3,iSkal)
        iBas_Eff = List_Bas(1,ilist_s)
        mdci  = iSD(10,iSkal)
        iShell= iSD(11,iSkal)
        index_i=list_bas(2,ilist_s)
        nFunc_i=iBas*iCmp

         Allocate(AO_vals(1:mAO,1:mGrid,1:iBas_Eff*iCmp))
         Call Get_AO_info(TabAO(ipTabAO(ilist_s)),mAO,mGrid,iBas_Eff,
     &         iCmp,AO_vals)

      Do iCB_Eff = 1, iBas_Eff*iCmp
         iCB = iCB_Eff + index_i - 1
!         iCB = Index(iCB_Eff)
         Do nGrad=1,3
           gf1(nGrad) = gf1(nGrad) +
     &                  CMO(iCB,i)*AO_vals(nGrad+1,iGrid,iCB_Eff)
         end do
      end do
      Do nGrad=1,3
        g_eff = list_g(nGrad,ilist_s)
        P2_ontop_d(g_eff,iGrid) = P2_ontop_d(g_eff,iGrid) +
     &    Fact*P2mo(ijkl)*gf1(nGrad)*TabMO(1,igrid,j)*
     &    TabMO(1,igrid,i)*TabMO(1,iGrid,l)
      end do
         Deallocate(AO_vals)
      end do

!Transform l to AO:
      gf1(1:3) = 0
      Do ilist_s=1,nlist_s
        iSkal=list_s(1,ilist_s)
        iCmp  = iSD( 2,iSkal)
        iBas  = iSD( 3,iSkal)
        iBas_Eff = List_Bas(1,ilist_s)
        mdci  = iSD(10,iSkal)
        iShell= iSD(11,iSkal)
        index_i=list_bas(2,ilist_s)
        nFunc_i=iBas*iCmp

         Allocate(AO_vals(1:mAO,1:mGrid,1:iBas_Eff*iCmp))
         Call Get_AO_info(TabAO(ipTabAO(ilist_s)),mAO,mGrid,iBas_Eff,
     &         iCmp,AO_vals)

      Do iCB_Eff = 1, iBas_Eff*iCmp
         iCB = iCB_Eff + index_i - 1
!         iCB = Index(iCB_Eff)
         Do nGrad=1,3
           gf1(nGrad) = gf1(nGrad) +
     &                  CMO(iCB,i)*AO_vals(nGrad+1,iGrid,iCB_Eff)
         end do
      end do
      Do nGrad=1,3
        g_eff = list_g(nGrad,ilist_s)
        P2_ontop_d(g_eff,iGrid) = P2_ontop_d(g_eff,iGrid) +
     &    Fact*P2mo(ijkl)*gf1(nGrad)*TabMO(1,igrid,j)*
     &    TabMO(1,igrid,k)*TabMO(1,iGrid,i)
      end do
         Deallocate(AO_vals)
      end do



*
                           End Do ! iGrid
                        End Do ! j_
                     End Do ! i_
                  End Do ! l_
               End Do ! k_
            End Do ! kIrrep
         End Do ! jIrrep
      End Do ! iIrrep
      End If
*
      RETURN
      END

         Subroutine Get_AO_info(TabAO,mAO,mGrid,iBas_eff,iCmp,AO_vals)
         Real*8 TabAO(mAO,mGrid,iBas_eff*iCmp)
         Real*8 AO_vals(mAO,mGrid,iBas_eff*iCmp)
         integer mAO,mGrid,iBas_Eff,iCmp

         AO_vals = TabAO
         End Subroutine
