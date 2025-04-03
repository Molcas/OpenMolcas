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

subroutine TimesE2_(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,ipCiOut)

use ipPage, only: ipin, opout, W
use MCLR_Data, only: n2Dens, nConf1, nDens2
use input_mclr, only: nRoots, nAsh, nRs2
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real*8 Kap(*)
integer ipCId, isym, jspin, ipS2, ipCiOut
real*8 ReCo
real*8 KapOut(*)
real*8 rdum(1)
real*8, allocatable :: RMOAA(:), Sc1(:), Sc2(:), Sc3(:), Temp4(:), Temp3(:)

call mma_allocate(RMOAA,n2dens,Label='RMOAA')
call mma_allocate(Sc1,ndens2,Label='Sc1')
call mma_allocate(Sc2,ndens2,Label='Sc2')
call mma_allocate(Sc3,ndens2,Label='Sc3')
call mma_allocate(Temp3,ndens2,Label='Temp3')
call mma_allocate(Temp4,ndens2,Label='Temp4')

if (doDMRG) then ! yma
  call dmrg_spc_change_mclr(RGras2(1:8),nash)
  call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
end if

call Uncompress(Kap,Sc1,isym)
call RInt_generic(SC1,rmoaa,rdum,Sc2,Temp3,Temp4,Sc3,isym,reco,jspin)

call Kap_CI(Temp4,nDens2,rmoaa,n2Dens,ipCIOUT)
call Ci_Ci(ipcid,ipS2)
call CI_KAP(ipCid,Sc1,Sc3,isym)

Sc1(:) = Sc2(:)+Sc3(:)

call Compress(Sc1,KapOut,isym)   ! ds
!call RecPrt('Ex',' ',KapOut,ndensC,1)

call ipin(ipS2)
call ipin(ipCIOUT)
W(ipCIOUT)%A(1:nConf1*nroots) = W(ipCIOUT)%A(1:nConf1*nroots)+W(ipS2)%A(1:nConf1*nroots)
call opOut(ipCId)

call mma_deallocate(Temp4)
call mma_deallocate(Temp3)
call mma_deallocate(Sc3)
call mma_deallocate(Sc2)
call mma_deallocate(Sc1)
call mma_deallocate(RMOAA)

if (doDMRG) call dmrg_spc_change_mclr(LRras2(1:8),nash)  ! yma

end subroutine TimesE2_
