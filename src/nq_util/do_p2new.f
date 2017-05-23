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
* Copyright (C) 2003, Sergey Gusarov                                   *
*               2003, Roland Lindh                                     *
*               2003, Per Ake Malmqvist                                *
************************************************************************
      Subroutine Do_P2new(P2mo,np2act,D1mo,nd1mo,TabMO,mAO,mGrid,nMOs,
     &                   P2_ontop,nP2_ontop,RhoI,RhoA,mRho)
************************************************************************
*                                                                      *
* Object: Calculation P2 ontop density and its derivatives             *
*                                                                      *
* Called from: Do_batch                                                *
*                                                                      *
* Calling    : FZero                                                   *
*                                                                      *
*     Authors: Sergey Gusarov, University of Lund, SWEDEN,             *
*              Roland Lindh,   University of Lund, SWEDEN,             *
*              P.-A.Malmqvist, University of Lund, SWEDEN,             *
*              Jan. 2003                                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "nq_info.fh"
#include "real.fh"
      Real*8 P2mo(np2act),D1mo(nd1mo),TabMO(mAO,mGrid,nMOs),
     &       P2_ontop(nP2_ontop,mGrid)
      Integer iOff_Ash(0:7), iOff_Bas(0:7), iOff_BasAct(0:7)
      Real*8 RhoI(mRho,mGrid)
      Real*8 RhoA(mRho,mGrid)
c
c      Common  /Coord/ GGrid(3,5000)
c      Common /iCoord/ mmGrid
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*      write(6,*) ' We are enter to P2cs'
*
*      Some paranoidal check(GS)
*
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
*                                                                      *
************************************************************************
*                                                                      *
*   P(1,...) - P_2                                                     *
*   P(2,...), P(3,...), P(4,...) - grad P_2                            *
*   P(5,...) - grad^2 P_2                                              *
*   P(6,...) - additional part grad^2 P_2 for CS functional            *
*                                                                      *
************************************************************************
*
*          Inactive part:
*
************************************************************************
*                                                                      *
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
c         Write(6,*)  " do_p2cs: Inact-Inact:", iIrrep,i,
c     &               TabMO(1,iGrid,i)
c
          RhoI(1,iGrid) = RhoI(1,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(1,iGrid,i)
          RhoI(2,iGrid) = RhoI(2,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(2,iGrid,i)
          RhoI(3,iGrid) = RhoI(3,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(3,iGrid,i)
          RhoI(4,iGrid) = RhoI(4,iGrid) +
     *                    TabMO(1,iGrid,i)*TabMO(4,iGrid,i)
*
          RhoI(5,iGrid) = RhoI(5,iGrid) +
     *                    TabMO( 2,iGrid,i)*TabMO( 2,iGrid,i) +
     *                    TabMO( 3,iGrid,i)*TabMO( 3,iGrid,i) +
     *                    TabMO( 4,iGrid,i)*TabMO( 4,iGrid,i)
*
          RhoI(6,iGrid) = RhoI(6,iGrid) + TabMO(1,iGrid,i)*
     *                    ( TabMO( 5,iGrid,i) +
     *                      TabMO( 8,iGrid,i) +
     *                      TabMO(10,iGrid,i) )
       End Do         ! i_
      End Do         ! iIrrep
      End Do         ! iGrid
*
      If (NumIsh.ne.0) Then
      Do iGrid = 1, mGrid
       P2_ontop(1,iGrid) = RhoI(1,iGrid)*RhoI(1,iGrid)
c      Write(6,'(4f28.20)') P2_ontop(1,iGrid)
*
       P2_ontop(2,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(2,iGrid)
*
       P2_ontop(3,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(3,iGrid)
*
       P2_ontop(4,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(4,iGrid)
*
       P2_ontop(5,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(5,iGrid) +
     *                     4.0d0*RhoI(1,iGrid)*RhoI(6,iGrid) +
     *                     8.0d0*(RhoI(2,iGrid)*RhoI(2,iGrid) +
     *                            RhoI(3,iGrid)*RhoI(3,iGrid) +
     *                            RhoI(4,iGrid)*RhoI(4,iGrid) )
*
       P2_ontop(6,iGrid) = 6.0d0*(RhoI(2,iGrid)*RhoI(2,iGrid) +
     *                            RhoI(3,iGrid)*RhoI(3,iGrid) +
     *                            RhoI(4,iGrid)*RhoI(4,iGrid))-
     *                      2.0d0*RhoI(1,iGrid)*RhoI(5,iGrid)
      End Do
      End If
************************************************************************
*
*          Active-Inactive part:
*
************************************************************************
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
                     RhoA(2,iGrid) = RhoA(2,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(2,iGrid,l)
                     RhoA(3,iGrid) = RhoA(3,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(3,iGrid,l)
                     RhoA(4,iGrid) = RhoA(4,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
                     RhoA(5,iGrid) = RhoA(5,iGrid) + D1mo(kl)*(
     *                   TabMO(2,iGrid,k)*TabMO(2,iGrid,l) +
     *                   TabMO(3,iGrid,k)*TabMO(3,iGrid,l) +
     *                   TabMO(4,iGrid,k)*TabMO(4,iGrid,l) )
                     RhoA(6,iGrid) = RhoA(6,iGrid) +
     *                   D1mo(kl)*TabMO(1,iGrid,k)*(
     *                   TabMO( 5,iGrid,k) +
     *                   TabMO( 8,iGrid,k) +
     *                   TabMO(10,iGrid,k) )
                  End Do     ! iGrid
               End Do      ! l_
            End Do       ! lIrrep
         End Do        ! k_
      End Do        ! kIrrep
*
      Do iGrid = 1, mGrid
       P2_ontop(1,iGrid) = P2_ontop(1,iGrid) +
     *                     RhoI(1,iGrid)*RhoA(1,iGrid)
       P2_ontop(2,iGrid) = P2_ontop(2,iGrid) +
     *                     2.0d0*RhoI(2,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(2,iGrid)
       P2_ontop(3,iGrid) = P2_ontop(3,iGrid) +
     *                     2.0d0*RhoI(3,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(3,iGrid)
       P2_ontop(4,iGrid) = P2_ontop(4,iGrid) +
     *                     2.0d0*RhoI(4,iGrid)*RhoA(1,iGrid) +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(4,iGrid)
       P2_ontop(5,iGrid) = P2_ontop(5,iGrid) +
     *                     2.0d0*RhoI(5,iGrid)*RhoA(1,iGrid)  +
     *                     2.0d0*RhoI(6,iGrid)*RhoA(1,iGrid)  +
     *                     8.0d0*(RhoI(2,iGrid)*RhoA(2,iGrid) +
     *                            RhoI(3,iGrid)*RhoA(3,iGrid) +
     *                            RhoI(4,iGrid)*RhoA(4,iGrid))+
     *                     2.0d0*RhoI(1,iGrid)*RhoA(5,iGrid)  +
     *                     2.0d0*RhoI(1,iGrid)*RhoA(6,iGrid)

*
       P2_ontop(6,iGrid) = P2_ontop(6,iGrid) +
     *                     6.0d0*(RhoI(2,iGrid)*RhoA(2,iGrid)  +
     *                            RhoI(3,iGrid)*RhoA(3,iGrid)  +
     *                            RhoI(4,iGrid)*RhoA(4,iGrid)) -
     *                            RhoI(5,iGrid)*RhoA(1,iGrid)  -
     *                            RhoI(1,iGrid)*RhoA(5,iGrid)
*
      End Do
*
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
*      Write(6,*) " do_p2cs:Irrep" ,iIrrep,jIrrep,kIrrep,ijkIrrep
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

                           Do iGrid = 1, mGrid
                              P2_ontop(1,iGrid) = P2_ontop(1,iGrid)  +
     *                                    Fact*P2mo(ijkl)*
     *                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     *                             TabMO(1,igrid,i)*TabMO(1,igrid,j)
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
                              Tmp1xx =          Fact*P2mo(ijkl)*(
     &                             TabMO(5,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(5,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(5,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(5,igrid,j)
     &                             )
                              Tmp1yy =          Fact*P2mo(ijkl)*(
     &                             TabMO(8,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(8,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(8,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(8,igrid,j)
     &                             )
                              Tmp1zz =          Fact*P2mo(ijkl)*(
     &                             TabMO(10,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(10,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(10,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(10,igrid,j)
     &                             )
c
                              Tmp2xx =          Fact*P2mo(ijkl)*(
     &                             TabMO(2,igrid,k)*TabMO(2,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(2,igrid,i)*TabMO(2,igrid,j)
     &                             )
                              Tmp2yy =          Fact*P2mo(ijkl)*(
     &                             TabMO(3,igrid,k)*TabMO(3,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(3,igrid,i)*TabMO(3,igrid,j)
     &                             )
                              Tmp2zz =          Fact*P2mo(ijkl)*(
     &                             TabMO(4,igrid,k)*TabMO(4,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(1,igrid,j) +
     &                             TabMO(1,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(4,igrid,i)*TabMO(4,igrid,j)
     &                             )
c
                              Tmp3xx =          Fact*P2mo(ijkl)*(
     &                             TabMO(2,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(2,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(2,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(2,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(2,igrid,l)*
     &                             TabMO(2,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(2,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(2,igrid,j)
     &                             )
                              Tmp3yy =          Fact*P2mo(ijkl)*(
     &                             TabMO(3,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(3,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(3,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(3,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(3,igrid,l)*
     &                             TabMO(3,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(3,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(3,igrid,j)
     &                             )
                              Tmp3zz =          Fact*P2mo(ijkl)*(
     &                             TabMO(4,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(4,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(4,igrid,k)*TabMO(1,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(4,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(4,igrid,l)*
     &                             TabMO(4,igrid,i)*TabMO(1,igrid,j) +
c
     &                             TabMO(1,igrid,k)*TabMO(4,igrid,l)*
     &                             TabMO(1,igrid,i)*TabMO(4,igrid,j)
     &                             )
                              P2_ontop(5,iGrid) = P2_ontop(5,iGrid) +
     &                                    Tmp1xx + Tmp1yy + Tmp1zz +
     &                             2.0d0*(Tmp2xx + Tmp2yy + Tmp2zz +
     &                                    Tmp3xx + Tmp3yy + Tmp3zz )
                              P2_ontop(6,iGrid) = P2_ontop(6,iGrid) +
     &                             Tmp3xx + Tmp3yy + Tmp3zz
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
