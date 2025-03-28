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

use ipPage, only: w
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use MCLR_Data, only: nConf1, n2Dens, nDens, nDens2
use input_mclr, only: nRoots, nAsh, nRS2
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2

implicit none
real*8 Kap(*)
integer ipCId, isym, jspin, ipS2, ipCiOut
real*8 reco
real*8 KapOut(*)
real*8 rdum(1)
real*8, allocatable :: Temp3(:), Temp4(:), Sc1(:), Sc2(:), Sc3(:), RMOAA(:)

call mma_allocate(RMOAA,n2Dens,Label='RMOAA')
call mma_allocate(Sc1,nDens2,Label='Sc1')
Sc1(:) = Zero
call mma_allocate(Sc2,nDens2,Label='Sc2')
call mma_allocate(Sc3,nDens2,Label='Sc3')
call mma_allocate(Temp3,nDens2,Label='Temp3')
call mma_allocate(Temp4,nDens2,Label='Temp4')

if (doDMRG) then ! yma
  call dmrg_spc_change_mclr(RGras2(1:8),nash)
  call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
end if
call Uncompress(Kap,Sc1,isym)

! Integral derivative !yma
call RInt_generic(SC1,rmoaa,rdum(1),Sc2,Temp3,Temp4,Sc3,isym,reco,jspin)

call Kap_CI(Temp4,nDens2,rmoaa,n2Dens,ipCIOUT)
call Ci_Ci(ipcid,ipS2)
call CI_KAP(ipCid,Sc1,Sc3,isym)

call DZaXpY(nDens,One,Sc2,1,Sc3,1,Sc1,1)

call Compress(Sc1,KapOut,isym)   ! ds
!call RecPrt('Ex',' ',KapOut,ndensC,1)

call ipin(ipS2)
call ipin(ipCIOUT)
call DaXpY_(nConf1*nroots,One,W(ipS2)%Vec,1,W(ipCIOUT)%Vec,1)
call opOut(ipCId)
! This is also orthogonalization of the solution vector
!do iR=1,nroots
!  do jR=1,nroots
!    ovl = ddot_(nconf1,W(ipciout)%Vec(iR),1,W(ipci)%Vec(jR),1)
!    call daxpy_(nconf1,-ovl,W(ipci)%Vec(jR),1,W(ipciout)%Vec(iR),1)
!  end do
!end do

call mma_deallocate(Temp4)
call mma_deallocate(Temp3)
call mma_deallocate(Sc3)
call mma_deallocate(Sc2)
call mma_deallocate(Sc1)
call mma_deallocate(rmoaa)

if (doDMRG) call dmrg_spc_change_mclr(LRras2(1:8),nash)  ! yma

end subroutine TimesE2
