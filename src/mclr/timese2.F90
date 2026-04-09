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

subroutine TimesE2(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,ipCiOut)

use ipPage, only: ipin, opout, W
use MCLR_Data, only: n2Dens, nConf1, nDens, nDensC
use input_mclr, only: nAsh, nRoots, nRS2
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Kap(nDensC), reco
integer(kind=iwp), intent(in) :: ipCId, isym, jspin, ipS2, ipCiOut
real(kind=wp), intent(out) :: KapOut(nDensC)
real(kind=wp) :: rdum(1)
real(kind=wp), allocatable :: RMOAA(:), Sc1(:), Sc2(:), Sc3(:), Temp3(:), Temp4(:)

call mma_allocate(RMOAA,n2Dens,Label='RMOAA')
call mma_allocate(Sc1,nDens,Label='Sc1')
Sc1(:) = Zero
call mma_allocate(Sc2,nDens,Label='Sc2')
call mma_allocate(Sc3,nDens,Label='Sc3')
call mma_allocate(Temp3,nDens,Label='Temp3')
call mma_allocate(Temp4,nDens,Label='Temp4')

if (doDMRG) then ! yma
  nash(:) = RGras2(:)
  nrs2(:) = RGras2(:)
end if
call Uncompress(Kap,Sc1,isym)

! Integral derivative !yma
call RInt_generic(SC1,rmoaa,rdum(1),Sc2,Temp3,Temp4,Sc3,isym,reco,jspin)

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
! This is also orthogonalization of the solution vector
!do iR=1,nroots
!  do jR=1,nroots
!    ovl = ddot_(nconf1,W(ipciout)%A(iR),1,W(ipci)%A(jR),1)
!    W(ipciout)%A(iR:iR+nconf1-1) = W(ipciout)%A(iR:iR+nconf1-1)-ovl*W(ipci)%A(jR:JR+nconf1-1)
!  end do
!end do

call mma_deallocate(Temp4)
call mma_deallocate(Temp3)
call mma_deallocate(Sc3)
call mma_deallocate(Sc2)
call mma_deallocate(Sc1)
call mma_deallocate(rmoaa)

if (doDMRG) nash(:) = LRras2(:)  ! yma

end subroutine TimesE2
