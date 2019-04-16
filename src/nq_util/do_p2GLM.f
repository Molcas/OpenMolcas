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
     &          RhoI,RhoA,mRho,nMOs,CMO,nAOs,nCMO,TabSO,nsym,ft)
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
#include "SysDef.fh"
#include "itmax.fh"
#include "info.fh"
#include "nq_info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
!#include "info.fh"
!Error could be TabAO...
      Integer iOff_Ash(0:7), iOff_Bas(0:7), iOff_BasAct(0:7),
     &        list_s(2,nlist_s),list_bas(2,nlist_s),Index(nIndex),
     &        list_g(3,nlist_s),ipTabAO(nlist_s,2),
     &        mAO,nAOs,mGrid,nP2_ontop,nGrad_Eff,np2act,nd1mo,nTabAO,
     &        mRho,nCMO,nsym
      Real*8 P2mo(np2act),D1mo(nd1mo),TabMO(mAO,mGrid,nMOs),
     &     P2_ontop(nP2_ontop,mGrid),TabAO(nTabAO),
     &     P2_ontop_d(np2_ontop,nGrad_Eff,mGrid),CMO(nCMO)
!      Real*8 P2AO(np2AO)
      logical ft
!      Real*8, allocatable, dimension(:,:,:) :: AO_vals
      Real*8, allocatable, dimension(:,:,:,:) :: dTabMO
!      Real*8, allocatable,dimension(:,:) :: CMO_u
!      Integer np2AO
      Real*8 RhoI(mRho,mGrid)
      Real*8 RhoA(mRho,mGrid)
      Real*8,dimension(1:mRho,1:mGrid,1:nGrad_Eff) :: dRhoI,dRhoA
!      Real*8 dRhoA(mRho,mGrid,nGrad_Eff)
!     Real*8 gf1(1:3)
      Logical Do_Grad
      integer g_eff,iGrid,CMO_OFF
      integer mo_i,ao_j
      Real*8 TabSO(mAO,mGrid,nMOs)
************************************************************************
*                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
      Call unused_logical(do_grad)
      Call unused_integer(naos)

!      write(*,*) 'CMO print'
!      do i=1,nMOs
!       write(*,*) 'MO number: ',i
!       do j=1,nAOs
!       write(*,*) j,CMO(j,i)
!       end do
!      end do

!unpack CMOs
!      nMO_IR = 0
!      nAO_IR = 0
!      Do iIrrep = 0, mIrrep-1
!        nMO_Ir = nMO_Ir + nIsh(iIrrep)+nAsh(iIrrep)
!        nAO_Ir = nAO_Ir + nBas(iIrrep)
!      End do
!      Allocate(CMO_u(1:nAO_IR,1:nMO_IR))
!        Call Fzero(CMO_u,nMO_IR*nAO_IR)



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
        iR = list_s(2,ilist_s)
        write(6,*) 'iR',iR
!        isym = NrOpr(iR,iOper,nSym)
!        write(6,*) 'isym',isym
        iAO = iSD(7,iSkal)
      end do
      end if



      If (nP2_ontop.eq.4) Then
         If (mAO.ne.10.or.mRho.ne.4) Then
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
      Call FZero(P2_ontop_d,nP2_ontop*nGrad_Eff*mGrid)
      dRhoI(1:mRho,1:mGrid,1:nGrad_Eff)=0.0d0
      dRhoA(1:mRho,1:mGrid,1:nGrad_Eff)=0.0d0
!      Call FZero(dRhoI,mRho*mGrid*nGrad_Eff)
!      Call FZero(dRhoA,mRho*mGrid*nGrad_Eff)
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

      !if(.true.) then
      !Assuming C1 symmetry right now.
      Allocate(dTabMO(1:nP2_ontop,1:nMOs,1:nGrad_eff,1:mgrid))
      dTabMO(1:nP2_ontop,1:nMOs,1:nGrad_eff,1:mGrid)=0.0d0

!      write(*,*) 'listbas 2'
!      do i=1,nlist_s
!      write(*,*) list_bas(2,i)
!      end do
!      write(*,*) 'nindex',nindex
!      do i=1,nindex
!      write(*,*) index(i)
!      end do
!      call xflush(6)
      ioff = 0

!      write(*,*) 'CMOs in p2glm'
!      do i=1,nCMO
!        write(*,*) CMO(i)
!      end do

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         Do ilist_s=1,nlist_s
            ish=list_s(1,ilist_s)
            iCmp  = iSD( 2,iSh)
            iBas  = iSD( 3,iSh)
            iBas_Eff = List_Bas(1,ilist_s)
            iPrim = iSD( 5,iSh)
            iAO   = iSD( 7,iSh)
            mdci  = iSD(10,iSh)
            iShell= iSD(11,iSh)

            kAO   = iCmp*iBas*mGrid
            nDeg  = nSym/nStab(mdci)
            nSO   = kAO*nDeg*mAO
            ipSOs = ipMem

            Call FZero(Work(ipSOs),nSO)

            iR=list_s(2,ilist_s)
            iSym=NrOpr(iR,iOper,nSym)

            Call SOAdpt_NQ(TabAO(ipTabAO(iList_s,1)),mAO,mGrid,iBas,
     &                  iBas_Eff,iCmp,iSym,Work(ipSOs),nDeg,iShell)

            Call GetMem('TmpCM','Allo','Real',ipTmpCMO,nCMO)
            Call GetMem('TDoIt','Allo','Inte',ipTDoIt,nMOs)
           !nPrint(135) = 99
           !do i=1,mAO*mGrid*iCmp*ndeg
           !end do
           Call FZero(TabSO,mAO*mGrid*nMOs)

            Call  SODist2(Work(ipSOs),mAO,mGrid,iBas,
     &                   iCmp,nDeg,TabSO,
     &                   iShell,nMOs,iAO,Work(ipTmpCMO),
     &                   nCMO,iWork(ipTDoIt))
            Call GetMem('TmpCM','Free','Real',ipTmpCMO,nCMO)
            Call GetMem('TDoIt','Free','Inte',ipTDoIt,nMOs)

            !Now I think we have the AOs converted to a full set of
            !SOs...
            MO_off = 0
            CMO_off = 0
            Do iIrrep = 0, mIrrep-1
              do MO_i = 1, nIsh(iIrrep)+nAsh(iIrrep)
                i = MO_i + MO_off
                do AO_j = 1, nBas(iIrrep)
                  j = AO_j + MO_off
                  do iCoord=1,3
                    g_eff = list_g(iCoord,ilist_s)
                    do iGrid = 1,mGrid

                  dTabMO(1,i,g_eff,iGrid) = dTabMO(1,i,g_eff,iGrid) +
     &            CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &            *TabSO(iCoord+1,iGrid,j)

!******************ADD STUFF FOR FTPBE HERE***************
        if (Functional_type.eq.GGA_type.and.ft) then

          Select CASE (iCoord)
           CASE (1) !x-coord
                    dTabMO(2,i,g_eff,iGrid) = dTabMO(2,i,g_eff,iGrid) +
     &              CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &             *(TabSO(5,iGrid,j))

                    dTabMO(3,i,g_eff,iGrid) = dTabMO(3,i,g_eff,iGrid) +
     &              CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &             *(TabSO(6,iGrid,j))

                    dTabMO(4,i,g_eff,iGrid) = dTabMO(4,i,g_eff,iGrid) +
     &              CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &            *(TabSO(7,iGrid,j))

           CASE (2) !y-coord
                    dTabMO(2,i,g_eff,iGrid) = dTabMO(2,i,g_eff,iGrid) +
     &              CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &             *(TabSO(6,iGrid,j))

                    dTabMO(3,i,g_eff,iGrid) = dTabMO(3,i,g_eff,iGrid) +
     &              CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &             *(TabSO(8,iGrid,j))

                    dTabMO(4,i,g_eff,iGrid) = dTabMO(4,i,g_eff,iGrid) +
     &              CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &            *(TabSO(9,iGrid,j))

           CASE (3) !z-coord
                    dTabMO(2,i,g_eff,iGrid) = dTabMO(2,i,g_eff,iGrid) +
     &              CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &             *(TabSO(7,iGrid,j))

                    dTabMO(3,i,g_eff,iGrid) = dTabMO(3,i,g_eff,iGrid) +
     &              CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &             *(TabSO(9,iGrid,j))

                    dTabMO(4,i,g_eff,iGrid) = dTabMO(4,i,g_eff,iGrid) +
     &              CMO(CMO_OFF + AO_j + (MO_i-1)*nBas(iIrrep))
     &            *(TabSO(10,iGrid,j))

           End Select

      end if
!*********************************************************

                    end do !grid
                  end do !coord
                end do!AO
              end do!MO
              MO_off = MO_off + nBas(iIrrep)
              CMO_off = CMO_off + nBas(iIrrep)**2
            end do!Irrep

      end do

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
          Do i_=1,nISh(iIrrep) + nFro(iIrrep)
            i = iOff_Bas(iIrrep) + i_
            RhoI(1,iGrid) = RhoI(1,iGrid) +
     &           TabMO(1,iGrid,i) * TabMO(1,iGrid,i)
        if (Functional_type.eq.GGA_type.and.ft) then
          RhoI(2,iGrid) = RhoI(2,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(2,iGrid,i)
          RhoI(3,iGrid) = RhoI(3,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(3,iGrid,i)
          RhoI(4,iGrid) = RhoI(4,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(4,iGrid,i)
        end if
*
      !Build dRhoI
        Do g_eff=1,nGrad_eff
              dRhoI(1,iGrid,g_eff) = dRhoI(1,iGrid,g_eff) +
     &        dTabMO(1,i,g_eff,iGrid)*TabMO(1,iGrid,i)!times 2 or not?

        if (Functional_type.eq.GGA_type.and.ft) then
            dRhoI(2,iGrid,g_eff) = dRhoI(2,iGrid,g_eff) +
     &      dTabMO(1,i,g_eff,iGrid)*TabMO(2,iGrid,i) +
     &      TabMO(1,iGrid,i)*dTabMO(2,i,g_eff,iGrid)

            dRhoI(3,iGrid,g_eff) = dRhoI(3,iGrid,g_eff) +
     &      dTabMO(1,i,g_eff,iGrid)*TabMO(3,iGrid,i) +
     &      TabMO(1,iGrid,i)*dTabMO(3,i,g_eff,iGrid)

            dRhoI(4,iGrid,g_eff) = dRhoI(4,iGrid,g_eff) +
     &      dTabMO(1,i,g_eff,iGrid)*TabMO(4,iGrid,i) +
     &      TabMO(1,iGrid,i)*dTabMO(4,i,g_eff,iGrid)

        end if !GGA
        end do !g_eff

          End Do         ! i_
        End Do         ! iIrrep
      End Do         ! iGrid
*
      If (NumIsh.ne.0) Then
      Do iGrid = 1, mGrid
         P2_ontop(1,iGrid) = RhoI(1,iGrid)*RhoI(1,iGrid)
*
        if (Functional_type.eq.GGA_type.and.ft) then
            P2_ontop(2,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(2,iGrid)
            P2_ontop(3,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(3,iGrid)
            P2_ontop(4,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(4,iGrid)
         end if
      do g_eff=1,nGrad_eff
        P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid) +
     & 4.0d0*dRhoI(1,iGrid,g_eff)*RhoI(1,iGrid)

        if (Functional_type.eq.GGA_type.and.ft) then
!******************ADD STUFF FOR FTPBE HERE***************
        P2_ontop_d(2,g_eff,iGrid) = P2_ontop_d(2,g_eff,iGrid) +
     &  4.0d0*dRhoI(2,iGrid,g_eff)*RhoI(1,iGrid) +
     & 8.0d0*dRhoI(1,iGrid,g_eff)*RhoI(2,iGrid)
!     &  4.0d0*dRhoI(1,iGrid,g_eff)*RhoI(2,iGrid) +
!     &  4.0d0*RhoI(1,iGrid)*dRhoI(2,iGrid,g_eff)

        P2_ontop_d(3,g_eff,iGrid) = P2_ontop_d(3,g_eff,iGrid) +
     &  4.0d0*dRhoI(3,iGrid,g_eff)*RhoI(1,iGrid) +
     & 8.0d0*dRhoI(1,iGrid,g_eff)*RhoI(3,iGrid)
!     &  4.0d0*dRhoI(1,iGrid,g_eff)*RhoI(3,iGrid) +
!     &  4.0d0*RhoI(1,iGrid)*dRhoI(3,iGrid,g_eff)

        P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid) +
     &  4.0d0*dRhoI(4,iGrid,g_eff)*RhoI(1,iGrid) +
     & 8.0d0*dRhoI(1,iGrid,g_eff)*RhoI(4,iGrid)
!     &  4.0d0*dRhoI(1,iGrid,g_eff)*RhoI(4,iGrid) +
!     &  4.0d0*RhoI(1,iGrid)*dRhoI(4,iGrid,g_eff)

        end if !GGA
      end do !ngrad


      End Do
      End If

************************************************************************
*          Active-Inactive part:                                       *
************************************************************************
      If (NumIsh.ne.0.and.NumAsh.ne.0) Then
       Do kIrrep = 0, mIrrep-1
         Do k_ = 1, nASh(kIrrep)
            k= k_ + iOff_BasAct(kIrrep)
            Do lIrrepx = 0, mIrrep-1
               Do l_ = 1, nAsh(lIrrepx)
                  l= l_ + iOff_BasAct(lIrrepx)
                  kl=iTri(k_ + iOff_Ash(kIrrep) ,
     &                    l_ + iOff_Ash(lIrrepx) )
                  Do iGrid = 1, mGrid
                     RhoA(1,iGrid) = RhoA(1,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(1,iGrid,l)
                     if (Functional_type.eq.GGA_type.and.ft) Then
                     RhoA(2,iGrid) = RhoA(2,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(2,iGrid,l)
                     RhoA(3,iGrid) = RhoA(3,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(3,iGrid,l)
                     RhoA(4,iGrid) = RhoA(4,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
                    end if

                    do g_eff=1,nGrad_eff
                      dRhoA(1,iGrid,g_eff) = dRhoA(1,iGrid,g_eff) +
     &                D1mo(kl)*dTabMO(1,k,g_eff,iGrid)*TabMO(1,iGrid,l)


!******************ADD STUFF FOR FTPBE HERE***************

                     if(Functional_type.eq.GGA_type.and.ft) Then
                       dRhoA(2,iGrid,g_eff) = dRhoA(2,iGrid,g_eff) +
     &               D1mo(kl)*dTabMO(1,k,g_eff,iGrid)*TabMO(2,iGrid,l) +
     &               D1mo(kl)*TabMO(1,iGrid,k)*dTabMO(2,l,g_eff,iGrid)

                       dRhoA(3,iGrid,g_eff) = dRhoA(3,iGrid,g_eff) +
     &               D1mo(kl)*dTabMO(1,k,g_eff,iGrid)*TabMO(3,iGrid,l) +
     &               D1mo(kl)*TabMO(1,iGrid,k)*dTabMO(3,l,g_eff,iGrid)

                       dRhoA(4,iGrid,g_eff) = dRhoA(4,iGrid,g_eff) +
     &               D1mo(kl)*dTabMO(1,k,g_eff,iGrid)*TabMO(4,iGrid,l) +
     &               D1mo(kl)*TabMO(1,iGrid,k)*dTabMO(4,l,g_eff,iGrid)
                     end if !GGA


                    end do
                  End Do     ! iGrid
               End Do      ! l_
            End Do       ! lIrrepx
         End Do        ! k_
       End Do        ! kIrrep
*
       Do iGrid = 1, mGrid
               P2_ontop(1,iGrid) = P2_ontop(1,iGrid) +
     *                           RhoI(1,iGrid)*RhoA(1,iGrid)
        if (Functional_type.eq.GGA_type.and.ft) Then
               P2_ontop(2,iGrid) = P2_ontop(2,iGrid) +
     *                     2.0d0*RhoI(2,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(2,iGrid)
               P2_ontop(3,iGrid) = P2_ontop(3,iGrid) +
     *                     2.0d0*RhoI(3,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(3,iGrid)
               P2_ontop(4,iGrid) = P2_ontop(4,iGrid) +
     *                     2.0d0*RhoI(4,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(4,iGrid)
        end if !gga
        do g_eff=1,nGrad_Eff
          P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid) +
     &                  2.0D0*RhoI(1,iGrid)*dRhoA(1,iGrid,g_eff) +
     &                  2.0D0*dRhoI(1,iGrid,g_eff)*RhoA(1,iGrid)
          if (Functional_type.eq.GGA_type.and.ft) Then

          P2_ontop_d(2,g_eff,iGrid) = P2_ontop_d(2,g_eff,iGrid) +
     &    2.0d0*dRhoI(2,iGrid,g_eff)*RhoA(1,iGrid) +
     &    4.0d0*RhoI(2,iGrid)*dRhoA(1,iGrid,g_eff) +
     &    4.0d0*dRhoI(1,iGrid,g_eff)*RhoA(2,iGrid) +
     &    2.0d0*RhoI(1,iGrid)*dRhoA(2,iGrid,g_eff)

          P2_ontop_d(3,g_eff,iGrid) = P2_ontop_d(3,g_eff,iGrid) +
     &    2.0d0*dRhoI(3,iGrid,g_eff)*RhoA(1,iGrid) +
     &    4.0d0*RhoI(3,iGrid)*dRhoA(1,iGrid,g_eff) +
     &    4.0d0*dRhoI(1,iGrid,g_eff)*RhoA(3,iGrid) +
     &    2.0d0*RhoI(1,iGrid)*dRhoA(3,iGrid,g_eff)
!     &    2.0d0*(dRhoI(3,iGrid,g_eff)*RhoA(1,iGrid) +
!     &           RhoI(3,iGrid)*dRhoA(1,iGrid,g_eff)) +
!     &    2.0d0*(dRhoI(1,iGrid,g_eff)*RhoA(3,iGrid) +
!     &           RhoI(1,iGrid)*dRhoA(3,iGrid,g_eff))

          P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid) +
     &    2.0d0*dRhoI(4,iGrid,g_eff)*RhoA(1,iGrid) +
     &    4.0d0*RhoI(4,iGrid)*dRhoA(1,iGrid,g_eff) +
     &    4.0d0*dRhoI(1,iGrid,g_eff)*RhoA(4,iGrid) +
     &    2.0d0*RhoI(1,iGrid)*dRhoA(4,iGrid,g_eff)
!     &    2.0d0*(dRhoI(4,iGrid,g_eff)*RhoA(1,iGrid) +
!     &           RhoI(4,iGrid)*dRhoA(1,iGrid,g_eff)) +
!     &    2.0d0*(dRhoI(1,iGrid,g_eff)*RhoA(4,iGrid) +
!     &           RhoI(1,iGrid)*dRhoA(4,iGrid,g_eff))
          end if !GGA
        end do !g_eff

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
                             if((ij.ge.kl).and.(k.eq.l)) Fact=1.0d0
                             if((kl.ge.ij).and.(i.eq.j)) Fact=1.0d0

                           Do iGrid = 1, mGrid
                             P2_ontop(1,iGrid) = P2_ontop(1,iGrid)  +
     *                                    Fact*P2mo(ijkl)*
     *                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     *                             TabMO(1,igrid,i)*TabMO(1,igrid,j)
                             if (Functional_type.eq.GGA_type.and.ft
     &                                     ) Then
                              P2_ontop(2,iGrid) = P2_ontop(2,iGrid)  +
     &                                     Fact*P2mo(ijkl)*(
     &                             TabMO(2,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
!c
     &                             TabMO(1,igrid,k)*TabMO(2,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
!c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(2,igrid,i)*TabMO(1,igrid,j) +
!c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(2,igrid,j)
     &                             )
                              P2_ontop(3,iGrid) = P2_ontop(3,iGrid) +
     &                                    Fact*P2mo(ijkl)*(
     &                             TabMO(3,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
!c
     &                             TabMO(1,igrid,k)*TabMO(3,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
!c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(3,igrid,i)*TabMO(1,igrid,j) +
!c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(3,igrid,j)
     &                             )
                              P2_ontop(4,iGrid) = P2_ontop(4,iGrid) +
     &                                    Fact*P2mo(ijkl)*(
     &                             TabMO(4,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
!c
     &                             TabMO(1,igrid,k)*TabMO(4,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
!c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(4,igrid,i)*TabMO(1,igrid,j) +
!c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(4,igrid,j)
     &                             )
                              end if
*
      do g_eff=1,nGrad_eff
        P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid) +
     &    Fact*P2mo(ijkl)*dTabMO(1,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &    TabMO(1,igrid,k)*TabMO(1,iGrid,l)

        P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid) +
     &    Fact*P2mo(ijkl)*dTabMO(1,j,g_eff,iGrid)*TabMO(1,igrid,i)*
     &    TabMO(1,igrid,k)*TabMO(1,iGrid,l)

        P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid) +
     &    Fact*P2mo(ijkl)*dTabMO(1,k,g_eff,iGrid)*TabMO(1,igrid,i)*
     &    TabMO(1,igrid,j)*TabMO(1,iGrid,l)

        P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid) +
     &    Fact*P2mo(ijkl)*dTabMO(1,l,g_eff,iGrid)*TabMO(1,igrid,i)*
     &    TabMO(1,igrid,j)*TabMO(1,iGrid,k)


         if (Functional_type.eq.GGA_type.and.ft) Then
           P2_ontop_d(2,g_eff,iGrid) = P2_ontop_d(2,g_eff,iGrid) +
     &     Fact*P2mo(ijkl)*
     &     (dTabMO(2,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &     TabMO(1,igrid,k)*TabMO(1,iGrid,l) +
     &     dTabMO(1,i,g_eff,iGrid)*TabMO(2,igrid,j)*
     &     TabMO(1,igrid,k)*TabMO(1,iGrid,l) +
     &     dTabMO(1,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &     TabMO(2,igrid,k)*TabMO(1,iGrid,l) +
     &     dTabMO(1,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &     TabMO(1,igrid,k)*TabMO(2,iGrid,l) +
     &
     &     TabMO(2,igrid,i)*dTabMO(1,j,g_eff,iGrid)*
     &     TabMO(1,iGrid,k)*TabMO(1,iGrid,l) +
     &     TabMO(1,igrid,i)*dTabMO(2,j,g_eff,iGrid)*
     &     TabMO(1,iGrid,k)*TabMO(1,iGrid,l) +
     &     TabMO(1,igrid,i)*dTabMO(1,j,g_eff,iGrid)*
     &     TabMO(2,iGrid,k)*TabMO(1,iGrid,l) +
     &     TabMO(1,igrid,i)*dTabMO(1,j,g_eff,iGrid)*
     &     TabMO(1,iGrid,k)*TabMO(2,iGrid,l) +
     &
     &     TabMO(2,iGrid,i)*TabMO(1,iGrid,j) *
     &     dTabMO(1,k,g_eff,iGrid)*TabMO(1,igrid,l) +
     &     TabMO(1,iGrid,i)*TabMO(2,iGrid,j) *
     &     dTabMO(1,k,g_eff,iGrid)*TabMO(1,igrid,l) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     dTabMO(2,k,g_eff,iGrid)*TabMO(1,igrid,l) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     dTabMO(1,k,g_eff,iGrid)*TabMO(2,igrid,l) +
     &
     &     TabMO(2,iGrid,i)*TabMO(1,iGrid,j) *
     &     TabMO(1,igrid,k)*dTabMO(1,l,g_eff,iGrid) +
     &     TabMO(1,iGrid,i)*TabMO(2,iGrid,j) *
     &     TabMO(1,igrid,k)*dTabMO(1,l,g_eff,iGrid) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     TabMO(2,igrid,k)*dTabMO(1,l,g_eff,iGrid) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     TabMO(1,igrid,k)*dTabMO(2,l,g_eff,iGrid))

           P2_ontop_d(3,g_eff,iGrid) = P2_ontop_d(3,g_eff,iGrid) +
     &     Fact*P2mo(ijkl)*
     &     (dTabMO(3,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &     TabMO(1,igrid,k)*TabMO(1,iGrid,l) +
     &     dTabMO(1,i,g_eff,iGrid)*TabMO(3,igrid,j)*
     &     TabMO(1,igrid,k)*TabMO(1,iGrid,l) +
     &     dTabMO(1,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &     TabMO(3,igrid,k)*TabMO(1,iGrid,l) +
     &     dTabMO(1,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &     TabMO(1,igrid,k)*TabMO(3,iGrid,l) +
     &
     &     TabMO(3,igrid,i)*dTabMO(1,j,g_eff,iGrid)*
     &     TabMO(1,iGrid,k)*TabMO(1,iGrid,l) +
     &     TabMO(1,igrid,i)*dTabMO(3,j,g_eff,iGrid)*
     &     TabMO(1,iGrid,k)*TabMO(1,iGrid,l) +
     &     TabMO(1,igrid,i)*dTabMO(1,j,g_eff,iGrid)*
     &     TabMO(3,iGrid,k)*TabMO(1,iGrid,l) +
     &     TabMO(1,igrid,i)*dTabMO(1,j,g_eff,iGrid)*
     &     TabMO(1,iGrid,k)*TabMO(3,iGrid,l) +
     &
     &     TabMO(3,iGrid,i)*TabMO(1,iGrid,j) *
     &     dTabMO(1,k,g_eff,iGrid)*TabMO(1,igrid,l) +
     &     TabMO(1,iGrid,i)*TabMO(3,iGrid,j) *
     &     dTabMO(1,k,g_eff,iGrid)*TabMO(1,igrid,l) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     dTabMO(3,k,g_eff,iGrid)*TabMO(1,igrid,l) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     dTabMO(1,k,g_eff,iGrid)*TabMO(3,igrid,l) +
     &
     &     TabMO(3,iGrid,i)*TabMO(1,iGrid,j) *
     &     TabMO(1,igrid,k)*dTabMO(1,l,g_eff,iGrid) +
     &     TabMO(1,iGrid,i)*TabMO(3,iGrid,j) *
     &     TabMO(1,igrid,k)*dTabMO(1,l,g_eff,iGrid) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     TabMO(3,igrid,k)*dTabMO(1,l,g_eff,iGrid) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     TabMO(1,igrid,k)*dTabMO(3,l,g_eff,iGrid))

           P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid) +
     &     Fact*P2mo(ijkl)*
     &     (dTabMO(4,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &     TabMO(1,igrid,k)*TabMO(1,iGrid,l) +
     &     dTabMO(1,i,g_eff,iGrid)*TabMO(4,igrid,j)*
     &     TabMO(1,igrid,k)*TabMO(1,iGrid,l) +
     &     dTabMO(1,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &     TabMO(4,igrid,k)*TabMO(1,iGrid,l) +
     &     dTabMO(1,i,g_eff,iGrid)*TabMO(1,igrid,j)*
     &     TabMO(1,igrid,k)*TabMO(4,iGrid,l) +
     &
     &     TabMO(4,igrid,i)*dTabMO(1,j,g_eff,iGrid)*
     &     TabMO(1,iGrid,k)*TabMO(1,iGrid,l) +
     &     TabMO(1,igrid,i)*dTabMO(4,j,g_eff,iGrid)*
     &     TabMO(1,iGrid,k)*TabMO(1,iGrid,l) +
     &     TabMO(1,igrid,i)*dTabMO(1,j,g_eff,iGrid)*
     &     TabMO(4,iGrid,k)*TabMO(1,iGrid,l) +
     &     TabMO(1,igrid,i)*dTabMO(1,j,g_eff,iGrid)*
     &     TabMO(1,iGrid,k)*TabMO(4,iGrid,l) +
     &
     &     TabMO(4,iGrid,i)*TabMO(1,iGrid,j) *
     &     dTabMO(1,k,g_eff,iGrid)*TabMO(1,igrid,l) +
     &     TabMO(1,iGrid,i)*TabMO(4,iGrid,j) *
     &     dTabMO(1,k,g_eff,iGrid)*TabMO(1,igrid,l) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     dTabMO(4,k,g_eff,iGrid)*TabMO(1,igrid,l) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     dTabMO(1,k,g_eff,iGrid)*TabMO(4,igrid,l) +
     &
     &     TabMO(4,iGrid,i)*TabMO(1,iGrid,j) *
     &     TabMO(1,igrid,k)*dTabMO(1,l,g_eff,iGrid) +
     &     TabMO(1,iGrid,i)*TabMO(4,iGrid,j) *
     &     TabMO(1,igrid,k)*dTabMO(1,l,g_eff,iGrid) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     TabMO(4,igrid,k)*dTabMO(1,l,g_eff,iGrid) +
     &     TabMO(1,iGrid,i)*TabMO(1,iGrid,j) *
     &     TabMO(1,igrid,k)*dTabMO(4,l,g_eff,iGrid))
         end if !GGA

      end do!loop over nGrad_eff
!  299 Continue
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
      !P2_ontop_d(1:nGrad_eff,1:mGrid) = 0


      !Do ilist_s=1,nlist_s !loop over AO shell lists
      !Do nGrad=1,3
      !  g_eff = list_g(nGrad,ilist_s)
      !end do
      !end do
      deAllocate(dTabMO)
      RETURN
      END subroutine

         Subroutine Get_AO_info(TabAO,mAO,mGrid,iBas_eff,iCmp,AO_vals)
         Real*8 TabAO(mAO,mGrid,iBas_eff*iCmp)
         Real*8 AO_vals(mAO,mGrid,iBas_eff*iCmp)
         integer mAO,mGrid,iBas_Eff,iCmp

         AO_vals = TabAO
         End Subroutine
