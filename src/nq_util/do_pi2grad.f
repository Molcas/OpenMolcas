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
      Subroutine Do_Pi2Grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &          P2_ontop,nP2_ontop,nGrad_Eff,
     &          list_s,nlist_s,list_bas,
     &          D1mo,nd1mo,TabMO,P2_ontop_d,
     &          RhoI,RhoA,mRho,nMOs,CMO,nCMO,TabSO,nsym,ft,
     &          P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,nPMO3p,MOs,
     &          MOx,MOy,MOz)
************************************************************************
*                                                                      *
* Object: Calculation P2 ontop density and its derivatives             *
*                                                                      *
* Called from: Do_batch                                                *
*                                                                      *
* Calling    : FZero                                                   *
*                                                                      *
*    INPUT:                                                            *
*   D1mo     = one-body density matrix in MO basis                     *
*   nd1mo    = size of D1mo                                            *
*   TabMO    = MO values computed on grid                              *
*   nMOs     = number of MO basis                                      *
*   mAO      = number of derivatives of AO...                          *
*   mGrid    = number of grid points                                   *
*                                                                      *
************************************************************************
      use iSD_data
      use Center_Info
      use Basis_Info, only: nBas
      use nq_pdft,    only: lft,lGGA
      use nq_Grid,    only: List_G
      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"
#include "nq_info.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
!Error could be TabAO...
      Integer list_s(2,nlist_s),list_bas(2,nlist_s),
     &        ipTabAO(nlist_s,2),
     &        mAO,nAOs,mGrid,nP2_ontop,nGrad_Eff,nd1mo,nTabAO,
     &        mRho,nCMO,nsym
      Real*8 D1mo(nd1mo),TabMO(mAO,mGrid,nMOs),
     &     P2_ontop(nP2_ontop,mGrid),TabAO(nTabAO),
     &     P2_ontop_d(np2_ontop,nGrad_Eff,mGrid),CMO(nCMO)
      logical ft
      Real*8, allocatable, dimension(:,:,:,:) :: dTabMO
      Real*8 RhoI(mRho,mGrid)
      Real*8 RhoA(mRho,mGrid)
      Real*8,dimension(1:mRho,1:mGrid,1:nGrad_Eff) :: dRhoI,dRhoA
      integer g_eff,iGrid
      Real*8 TabSO(mAO,mGrid,nMOs)
      Real*8,DIMENSION(mGrid*NASHT)::P2MOCube,MOs,dMOs,MOx,MOy,MOz
      INTEGER IOff1,iOff2,iOff3,nPi,iCoordOff,iGridOff,iCoord,nBasf,
     &        nOccO,nPMO3p,
     &        iCoord1,iCoord2,iCoord3,iCoordOff1,iCoordOff2,iCoordOff3
      Real*8,DIMENSION(nPMO3p)::P2MOCubex,P2MOCubey,P2MOCubez,
     &                          dMOx,dMOy,dMOz

      Real*8 TabSO2(mAO*mGrid*nMOs)
      Real*8 dTabMO2(nMOs)

************************************************************************
*                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
      Call unused_integer(naos)

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
*   P(1,...) - P_2                                                     *
*   P(2,...), P(3,...), P(4,...) - grad P_2                            *
************************************************************************

      Call mma_Allocate(dTabMO,nP2_ontop,nMOs,nGrad_eff,mgrid,
     &                  Label='dTabMO')
      dTabMO(:,:,:,:)=Zero

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Do ilist_s=1,nlist_s
         ish=list_s(1,ilist_s)
         iCmp  = iSD( 2,iSh)
         iBas  = iSD( 3,iSh)
         iBas_Eff = List_Bas(1,ilist_s)
         iAO   = iSD( 7,iSh)
         mdci  = iSD(10,iSh)

         nDeg  = nSym/dc(mdci)%nStab

         Call FZero(TabSO,mAO*mGrid*nMOs)

         iR=list_s(2,ilist_s)
         iSym=NrOpr(iR)

         Call SOAdpt_NQ(TabAO(ipTabAO(iList_s,1)),mAO,mGrid,iBas,
     &                  iBas_Eff,iCmp,iSym,TabSO,nMOs,nDeg,iAO)

         CALL ConvertTabSO(TabSO2,TabSO,mAO,mGrid,nMOs)

         Do iGrid=1,mGrid
            IGridOff=(iGrid-1)*mAO*nMOs
            Do iCoord=1,3
               ICoordOff=IGridOff+(iCoord-1)*nMOs
               g_eff = list_g(iCoord,ilist_s)
               ICoordOff1=0
               ICoordOff2=0
               ICoordOff3=0
               IF (lft.and.lGGA) THEN
                  If (iCoord.eq.1) Then
                     iCoord1=4
                     iCoord2=5
                     iCoord3=6
                  Else If(iCoord.eq.2) Then
                     iCoord1=5
                     iCoord2=7
                     iCoord3=8
                  Else If(iCoord.eq.3) Then
                     iCoord1=6
                     iCoord2=8
                     iCoord3=9
                  End If
                  ICoordOff1=IGridOff+(iCoord1-1)*nMOs
                  ICoordOff2=IGridOff+(iCoord2-1)*nMOs
                  ICoordOff3=IGridOff+(iCoord3-1)*nMOs
               END IF
               do iIrrep=0,mIrrep-1
                  nOccO=nIsh(iIrrep)+nAsh(iIrrep)
                  IF (nOccO.eq.0) CYCLE
                  nBasF=nBas(iIrrep)

                  CALL DGEMM_('T','N',nOccO,1,nBasF,
     &                     1.0d0,CMO(OffBas2(iIrrep)),nBasF,
     &                           TabSO2(iCoordOff+OffBas(iIrrep)),nBasF,
     &                     0.0d0,dTabMO2,nOccO)
                  CALL DAXPY_(nOccO,1.0d0,dTabMO2,1,
     &      dTabMO(1,OffBas(iIrrep)+nFro(iIrrep),g_eff,iGrid),nP2_ontop)

                  IF (lft.and.lGGA) THEN
                     CALL DGEMM_('T','N',nOccO,1,nBasF,
     &                    1.0d0,CMO(OffBas2(iIrrep)),nBasF,
     &                          TabSO2(iCoordOff1+OffBas(iIrrep)),nBasF,
     &                    0.0d0,dTabMO2,nOccO)
                     CALL DAXPY_(nOccO,1.0d0,dTabMO2,1,
     &      dTabMO(2,OffBas(iIrrep)+nFro(iIrrep),g_eff,iGrid),nP2_ontop)

                     CALL DGEMM_('T','N',nOccO,1,nBasF,
     &                    1.0d0,CMO(OffBas2(iIrrep)),nBasF,
     &                          TabSO2(iCoordOff2+OffBas(iIrrep)),nBasF,
     &                    0.0d0,dTabMO2,nOccO)
                     CALL DAXPY_(nOccO,1.0d0,dTabMO2,1,
     &      dTabMO(3,OffBas(iIrrep)+nFro(iIrrep),g_eff,iGrid),nP2_ontop)

                     CALL DGEMM_('T','N',nOccO,1,nBasF,
     &                    1.0d0,CMO(OffBas2(iIrrep)),nBasF,
     &                          TabSO2(iCoordOff3+OffBas(iIrrep)),nBasF,
     &                    0.0d0,dTabMO2,nOccO)
                     CALL DAXPY_(nOccO,1.0d0,dTabMO2,1,
     &      dTabMO(4,OffBas(iIrrep)+nFro(iIrrep),g_eff,iGrid),nP2_ontop)
                  END IF
               end do ! iIrrep
            End Do  ! iCoord
         End Do   ! iGrid
      END DO         ! iList_s
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
!******************ADD STUFF FOR FT: HERE***************
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


!******************ADD STUFF FOR FT: HERE***************

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

          P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid) +
     &    2.0d0*dRhoI(4,iGrid,g_eff)*RhoA(1,iGrid) +
     &    4.0d0*RhoI(4,iGrid)*dRhoA(1,iGrid,g_eff) +
     &    4.0d0*dRhoI(1,iGrid,g_eff)*RhoA(4,iGrid) +
     &    2.0d0*RhoI(1,iGrid)*dRhoA(4,iGrid,g_eff)
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
       nPi=nP2_ontop
*
       DO g_eff=1,nGrad_eff
        Do iGrid=1,mGrid
         IOff1=(iGrid-1)*NASHT
         do iIrrep=0,mIrrep-1
          IOff2=IOff_Ash(iIrrep)+1
          IOff3=IOff_BasAct(iIrrep)+1
          CALL DCopy_(nAsh(iIrrep),dTabMO(1,iOff3,g_eff,iGrid),nPi,
     &                               dMOs(IOff1+IOff2)        ,1  )
          if(lft.and.lGGA) then
           CALL DCopy_(nAsh(iIrrep),dTabMO(2,iOff3,g_eff,iGrid),nPi,
     &                                dMOx(IOff1+IOff2)        ,1  )
           CALL DCopy_(nAsh(iIrrep),dTabMO(3,iOff3,g_eff,iGrid),nPi,
     &                                dMOy(IOff1+IOff2)        ,1  )
           CALL DCopy_(nAsh(iIrrep),dTabMO(4,iOff3,g_eff,iGrid),nPi,
     &                                dMOz(IOff1+IOff2)        ,1  )
          end if
         end do
        End Do

        Do iGrid=1,mGrid
         IOff1=(iGrid-1)*NASHT
         do IIrrep=0,mIrrep-1
          IOff2=IOff1+iOff_Ash(IIrrep)+1
          P2_ontop_d(1,g_eff,iGrid)=P2_ontop_d(1,g_eff,iGrid)+
     &    4.0d0*ddot_(nAsh(IIrrep),dMOs(IOff2),1,P2MOCube(IOff2),1)
          if(lft.and.lGGA) then
           P2_ontop_d(2,g_eff,iGrid)=P2_ontop_d(2,g_eff,iGrid)+
     &     4.0d0*ddot_(nAsh(IIrrep),dMOx(IOff2),1,P2MOCube(IOff2),1)+
     &     4.0d0*ddot_(nAsh(IIrrep),dMOs(IOff2),1,P2MOCubex(IOff2),1)
           P2_ontop_d(3,g_eff,iGrid)=P2_ontop_d(3,g_eff,iGrid)+
     &     4.0d0*ddot_(nAsh(IIrrep),dMOy(IOff2),1,P2MOCube(IOff2),1)+
     &     4.0d0*ddot_(nAsh(IIrrep),dMOs(IOff2),1,P2MOCubey(IOff2),1)
           P2_ontop_d(4,g_eff,iGrid)=P2_ontop_d(4,g_eff,iGrid)+
     &     4.0d0*ddot_(nAsh(IIrrep),dMOz(IOff2),1,P2MOCube(IOff2),1)+
     &     4.0d0*ddot_(nAsh(IIrrep),dMOs(IOff2),1,P2MOCubez(IOff2),1)
          end if
         end do
        End Do
       END DO


       DO iGrid=1,mGrid
        IOff1=(iGrid-1)*NASHT
        Do IIrrep=0,mIrrep-1
         IOff2=IOff1+iOff_Ash(IIrrep)+1
         P2_ontop(1,iGrid)=P2_ontop(1,iGrid)+
     &   ddot_(nAsh(IIrrep),MOs(IOff2),1,P2MOCube(IOff2),1)
        End Do
       END DO

       If(lGGA.and.lft) Then
        DO iGrid=1,mGrid
         IOff1=(iGrid-1)*NASHT
         Do IIrrep=0,mIrrep-1
          IOff2=IOff1+iOff_Ash(IIrrep)+1
          P2_ontop(2,iGrid)=P2_ontop(2,iGrid)+
     &    4.0d0*ddot_(nAsh(IIrrep),MOx(IOff2),1,P2MOCube(IOff2),1)
          P2_ontop(3,iGrid)=P2_ontop(3,iGrid)+
     &    4.0d0*ddot_(nAsh(IIrrep),MOy(IOff2),1,P2MOCube(IOff2),1)
          P2_ontop(4,iGrid)=P2_ontop(4,iGrid)+
     &    4.0d0*ddot_(nAsh(IIrrep),MOz(IOff2),1,P2MOCube(IOff2),1)
C          write(6,*)'MOz used for dPiz'
C          write(6,*) iGrid,iOff2-iOff1
C          CALL RecPrt(' ',' ',MOz(iOff2),1,nAsh(iIrrep))
         End Do
        END DO
       End If

      End If
      Call mma_deAllocate(dTabMO)
      RETURN
      END subroutine

