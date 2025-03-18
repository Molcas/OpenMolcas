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
! Copyright (C) 2017, Andrew M. Sand                                   *
!                                                                      *
!***********************************************************************
      Subroutine TimesE2_(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,ipCiOut)
!
      use ipPage, only: W
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: One
      use MCLR_Data
      use input_mclr, only: nRoots,nAsh,nRs2
      use dmrginfo, only: DoDMRG,LRRAS2,RGRAS2
      Implicit None
      Real*8 Kap(*)
      Integer ipCId,isym,jspin,ipS2,ipCiOut
      Real*8 ReCo
      Real*8 KapOut(*)

      Integer opOut
      Real*8 rdum(1)
      Real*8, Allocatable:: RMOAA(:), Sc1(:), Sc2(:), Sc3(:),           &
     &                      Temp4(:), Temp3(:)
      Integer iRC
      Integer, External:: ipIN
!
      Call mma_allocate(RMOAA,n2dens,Label='RMOAA')
      Call mma_allocate(Sc1,ndens2,Label='Sc1')
      Call mma_allocate(Sc2,ndens2,Label='Sc2')
      Call mma_allocate(Sc3,ndens2,Label='Sc3')
      Call mma_allocate(Temp3,ndens2,Label='Temp3')
      Call mma_allocate(Temp4,ndens2,Label='Temp4')
!
      if(doDMRG)then ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if

      Call Uncompress(Kap,Sc1,isym)
      Call RInt_generic(SC1,rmoaa,rdum,Sc2,Temp3,Temp4,Sc3,             &
     &                 isym,reco,jspin)

      Call Kap_CI(Temp4,nDens2,rmoaa,n2Dens,ipCIOUT)
      Call Ci_Ci(ipcid,ipS2)
      Call CI_KAP(ipCid,Sc1,Sc3,isym)

      Call DZaXpY(nDens,One,Sc2,1,Sc3,1,Sc1,1)
!
      Call Compress(Sc1,KapOut,isym)   ! ds
!     Call RecPrt('Ex',' ',KapOut,ndensC,1)
!
      irc=ipin(ipS2)
      irc=ipin(ipCIOUT)
      Call DaXpY_(nConf1*nroots,One,W(ipS2)%Vec,1,W(ipCIOUT)%Vec,1)
      irc=opOut(ipCId)

!
      Call mma_deallocate(Temp4)
      Call mma_deallocate(Temp3)
      Call mma_deallocate(Sc3)
      Call mma_deallocate(Sc2)
      Call mma_deallocate(Sc1)
      Call mma_deallocate(RMOAA)

      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
      end if
!
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(irc)
#endif
      End Subroutine TimesE2_
