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
      Subroutine Do_RhoIA(nRho,mGrid,Rho,RhoI,
     &                     RhoA,mRho,TabMO,mAO,nMOs,
     &                     D1mo,nd1mo)
      Implicit Real*8 (a-h,o-z)
      Dimension RhoI(mRho,mGrid),RhoA(mRho,mGrid)
      Dimension TabMO(mAO,mGrid,nMOs)
      Dimension Rho(nRho,mGrid)
      Dimension D1mo(nd1mo)
      Integer iOff_Ash(0:7), iOff_Bas(0:7), iOff_BasAct(0:7)
#include "nq_info.fh"
#include "real.fh"
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      Call FZero(RhoI,mRho*mGrid)
      Call FZero(RhoA,mRho*mGrid)
      jOffA_ = 0
      jOffB_ = 0
      Do iIrrep = 0, mIrrep-1
         iOff_Ash(iIrrep)=jOffA_
         iOff_Bas(iIrrep)=jOffB_
         iOff_BasAct(iIrrep)=jOffB_ + nIsh(iIrrep) + nFro(iIrrep)
         jOffA_=jOffA_+nAsh(iIrrep)
         jOffB_=jOffB_+mBas(iIrrep)
      End Do
*
      NumIsh = 0
      NumAsh = 0
      Do iIrrep=0, mIrrep-1
         NumIsh = NumIsh + nISh(iIrrep)
         NumAsh = NumAsh + nAsh(iIrrep)
      End Do

      If (mRho.eq.1) Then
      Do iGrid = 1, mGrid
      Do iIrrep=0, mIrrep-1
*       Write(6,*) ' nISh',nISh(iIrrep)
       Do i_=1,nISh(iIrrep) + nFro(iIrrep)
        i = iOff_Bas(iIrrep) + i_
          RhoI(1,iGrid) = RhoI(1,iGrid) +
     &                    TabMO(1,iGrid,i)*TabMO(1,iGrid,i)
       End Do         ! i_
      End Do         ! iIrrep
      End Do         ! iGrid
*
      Do kIrrep = 0, mIrrep-1
         Do k_ = 1, nASh(kIrrep)
            k= k_ + iOff_BasAct(kIrrep)
            Do lIrrep = 0, mIrrep-1
               Do l_ = 1, nAsh(lIrrep)
                  l= l_ + iOff_BasAct(lIrrep)
                  kl=iTri(k_ + iOff_Ash(kIrrep) ,
     &                    l_ + iOff_Ash(lIrrep) )
                  Do iGrid = 1, mGrid
                     RhoA(1,iGrid) = RhoA(1,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(1,iGrid,l)
                  End Do     ! iGrid
               End Do      ! l_
            End Do       ! lIrrep
         End Do        ! k_
      End Do        ! kIrrep
      Else If (mRho.eq.4) Then
      Do iGrid = 1, mGrid
      Do iIrrep=0, mIrrep-1
*       Write(6,*) ' nISh',nISh(iIrrep)
       Do i_=1,nISh(iIrrep) + nFro(iIrrep)
        i = iOff_Bas(iIrrep) + i_
          RhoI(1,iGrid) = RhoI(1,iGrid) +
     &                    TabMO(1,iGrid,i)*TabMO(1,iGrid,i)
          RhoI(2,iGrid) = RhoI(2,iGrid) +
     &                    TabMO(1,iGrid,i)*TabMO(2,iGrid,i)
          RhoI(3,iGrid) = RhoI(3,iGrid) +
     &                    TabMO(1,iGrid,i)*TabMO(3,iGrid,i)
          RhoI(4,iGrid) = RhoI(4,iGrid) +
     &                    TabMO(1,iGrid,i)*TabMO(4,iGrid,i)
       End Do         ! i_
      End Do         ! iIrrep
      End Do         ! iGrid
*
      Do kIrrep = 0, mIrrep-1
         Do k_ = 1, nASh(kIrrep)
            k= k_ + iOff_BasAct(kIrrep)
            Do lIrrep = 0, mIrrep-1
               Do l_ = 1, nAsh(lIrrep)
                  l= l_ + iOff_BasAct(lIrrep)
                  kl=iTri(k_ + iOff_Ash(kIrrep) ,
     &                    l_ + iOff_Ash(lIrrep) )
                  Do iGrid = 1, mGrid
                     RhoA(1,iGrid) = RhoA(1,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(1,iGrid,l)
                     RhoA(2,iGrid) = RhoA(2,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(2,iGrid,l)
                     RhoA(3,iGrid) = RhoA(3,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(3,iGrid,l)
                     RhoA(4,iGrid) = RhoA(4,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
                  End Do     ! iGrid
               End Do      ! l_
            End Do       ! lIrrep
         End Do        ! k_
      End Do        ! kIrrep
      Else If (mRho.eq.6) Then
      Do iGrid = 1, mGrid
      Do iIrrep=0, mIrrep-1
       Do i_=1,nISh(iIrrep) + nFro(iIrrep)
        i = iOff_Bas(iIrrep) + i_
          RhoI(1,iGrid) = RhoI(1,iGrid) +
     &                    TabMO(1,iGrid,i)*TabMO(1,iGrid,i)
          RhoI(2,iGrid) = RhoI(2,iGrid) +
     &                    TabMO(1,iGrid,i)*TabMO(2,iGrid,i)
          RhoI(3,iGrid) = RhoI(3,iGrid) +
     &                    TabMO(1,iGrid,i)*TabMO(3,iGrid,i)
          RhoI(4,iGrid) = RhoI(4,iGrid) +
     &                    TabMO(1,iGrid,i)*TabMO(4,iGrid,i)
*
          RhoI(5,iGrid) = RhoI(5,iGrid) +
     &                    TabMO( 2,iGrid,i)*TabMO( 2,iGrid,i) +
     &                    TabMO( 3,iGrid,i)*TabMO( 3,iGrid,i) +
     &                    TabMO( 4,iGrid,i)*TabMO( 4,iGrid,i)
*
          RhoI(6,iGrid) = RhoI(6,iGrid) + TabMO(1,iGrid,i)*
     &                    ( TabMO( 5,iGrid,i) +
     &                      TabMO( 8,iGrid,i) +
     &                      TabMO(10,iGrid,i) )
       End Do         ! i_
      End Do         ! iIrrep
      End Do         ! iGrid
*
      Do kIrrep = 0, mIrrep-1
         Do k_ = 1, nASh(kIrrep)
            k= k_ + iOff_BasAct(kIrrep)
            Do lIrrep = 0, mIrrep-1
               Do l_ = 1, nAsh(lIrrep)
                  l= l_ + iOff_BasAct(lIrrep)
                  kl=iTri(k_ + iOff_Ash(kIrrep) ,
     &                    l_ + iOff_Ash(lIrrep) )
                  Do iGrid = 1, mGrid
                     RhoA(1,iGrid) = RhoA(1,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(1,iGrid,l)
                     RhoA(2,iGrid) = RhoA(2,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(2,iGrid,l)
                     RhoA(3,iGrid) = RhoA(3,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(3,iGrid,l)
                     RhoA(4,iGrid) = RhoA(4,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
                     RhoA(5,iGrid) = RhoA(5,iGrid) + D1mo(kl)*(
     &                   TabMO(2,iGrid,k)*TabMO(2,iGrid,l) +
     &                   TabMO(3,iGrid,k)*TabMO(3,iGrid,l) +
     &                   TabMO(4,iGrid,k)*TabMO(4,iGrid,l) )
                     RhoA(6,iGrid) = RhoA(6,iGrid) +
     &                   D1mo(kl)*TabMO(1,iGrid,k)*(
     &                   TabMO( 5,iGrid,k) +
     &                   TabMO( 8,iGrid,k) +
     &                   TabMO(10,iGrid,k) )
                  End Do     ! iGrid
               End Do      ! l_
            End Do       ! lIrrep
         End Do        ! k_
      End Do        ! kIrrep
      Else
         Call WarningMessage(2,
     &              ' Somethings wrong in RhoI,A calculation')
         Call Abend()
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Rho)
      End
