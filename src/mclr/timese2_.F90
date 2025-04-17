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
use MCLR_Data, only: n2Dens, nConf1, nDens, nDensC
use input_mclr, only: nAsh, nRoots, nRs2
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Kap(nDensC), ReCo
integer(kind=iwp), intent(in) :: ipCId, isym, jspin, ipS2, ipCiOut
real(kind=wp), intent(out) :: KapOut(nDensC)
real(kind=wp) :: rdum(1)
real(kind=wp), allocatable :: RMOAA(:), Sc1(:), Sc2(:), Sc3(:), Temp3(:), Temp4(:)

call mma_allocate(RMOAA,n2Dens,Label='RMOAA')
call mma_allocate(Sc1,nDens,Label='Sc1')
call mma_allocate(Sc2,nDens,Label='Sc2')
call mma_allocate(Sc3,nDens,Label='Sc3')
call mma_allocate(Temp3,nDens,Label='Temp3')
call mma_allocate(Temp4,nDens,Label='Temp4')

if (doDMRG) then ! yma
  nash(:) = RGras2(:)
  nrs2(:) = RGras2(:)
end if

call Uncompress(Kap,Sc1,isym)
call RInt_generic(SC1,rmoaa,rdum,Sc2,Temp3,Temp4,Sc3,isym,reco,jspin)

call Kap_CI(Temp4,nDens,rmoaa,n2Dens,ipCIOUT)
call Ci_Ci(ipcid,ipS2)
call CI_KAP(ipCid,Sc1,Sc3,isym)

Sc1(:) = Sc2(:)+Sc3(:)

call Compress(Sc1,KapOut,isym)   ! ds
!call RecPrt('Ex',' ',KapOut,nDensC,1)

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

if (doDMRG) nash(:) = LRras2(:) ! yma

end subroutine TimesE2_
