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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************
      Subroutine TimesE2(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,ipCiOut)
      use ipPage, only: w
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, One
      use MCLR_Data, only: nConf1,n2Dens,nDens,nDens2
      use input_mclr, only: nRoots,nAsh,nRS2
      use dmrginfo, only: DoDMRG,LRRAS2,RGRAS2
      Implicit None
      Real*8 Kap(*)
      Integer ipCId,isym,jspin,ipS2,ipCiOut
      Real*8 reco
      Real*8 KapOut(*)

      Integer opOut
      Real*8 rdum(1)
      Real*8, Allocatable:: Temp3(:), Temp4(:),                         &
     &                      Sc1(:), Sc2(:), Sc3(:), RMOAA(:)
      Integer iRC
      Integer, External:: ipIn
!
      Call mma_allocate(RMOAA,n2Dens,Label='RMOAA')
      Call mma_allocate(Sc1,nDens2,Label='Sc1')
      Sc1(:)=Zero
      Call mma_allocate(Sc2,nDens2,Label='Sc2')
      Call mma_allocate(Sc3,nDens2,Label='Sc3')
      Call mma_allocate(Temp3,nDens2,Label='Temp3')
      Call mma_allocate(Temp4,nDens2,Label='Temp4')
!
      if(doDMRG)then ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if
      Call Uncompress(Kap,Sc1,isym)

! Integral derivative !yma
      Call RInt_generic(SC1,rmoaa,rdum(1),                              &
     &                 Sc2,                                             &
     &                 Temp3,Temp4,Sc3,                                 &
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
      Call DaXpY_(nConf1*nroots,One,                                    &
     &               W(ipS2)%Vec,1,                                     &
     &               W(ipCIOUT)%Vec,1)
      irc=opOut(ipCId)
      !! This is also orthogonalization of the solution vector
!     do iR = 1, nroots
!       do jR = 1, nroots
!         ovl = ddot_(nconf1,work(ipin(ipciout)+(iR-1)*nconf1),1,
!    *                       work(ipin(ipci)+(jR-1)*nconf1),1)
!         call daxpy_(nconf1,-ovl,work(ipin(ipci)+(jR-1)*nconf1),1,
!    *                            work(ipin(ipciout)+(iR-1)*nconf1),1)
!       end do
!     end do

!
      Call mma_deallocate(Temp4)
      Call mma_deallocate(Temp3)
      Call mma_deallocate(Sc3)
      Call mma_deallocate(Sc2)
      Call mma_deallocate(Sc1)
      Call mma_deallocate(rmoaa)

      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
      end if
!
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(irc)
#endif
      End Subroutine TimesE2
