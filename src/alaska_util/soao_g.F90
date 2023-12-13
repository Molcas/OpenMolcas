!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine SOAO_g(iSD4,nSD,nSO,MemPrm,MemMax,iBsInc,jBsInc,kBsInc,lBsInc,iPrInc,jPrInc,kPrInc,lPrInc,ipMem1,ipMem2,Mem1,Mem2,iFnc, &
                  MemPSO)

use Basis_Info, only: Shells
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4), nSO, MemPrm, MemMax, ipMem1
integer(kind=iwp), intent(out) :: iBsInc, jBsInc, kBsInc, lBsInc, iPrInc, jPrInc, kPrInc, lPrInc, ipMem2, Mem1, Mem2, iFnc(4), &
                                  MemPSO
integer(kind=iwp) :: iAnga(4), iAO(4), iBasi, iCmpa(4), iPrimi, iQuad, jBasj, jPrimj, kBask, kPrimk, lBasl, lPriml

iPrimi = Shells(iSD4(0,1))%nExp
jPrimj = Shells(iSD4(0,2))%nExp
kPrimk = Shells(iSD4(0,3))%nExp
lPriml = Shells(iSD4(0,4))%nExp
iBasi = Shells(iSD4(0,1))%nBasis
jBasj = Shells(iSD4(0,2))%nBasis
kBask = Shells(iSD4(0,3))%nBasis
lBasl = Shells(iSD4(0,4))%nBasis
do iQuad=1,4
  iAnga(iQuad) = iSD4(1,iQuad)
  iCmpa(iQuad) = iSD4(2,iQuad)
  iAO(iQuad) = iSD4(7,iQuad)
end do

call PSOAO1(nSO,MemPrm,MemMax,iAnga,iCmpa,iAO,iFnc,iBasi,iBsInc,jBasj,jBsInc,kBask,kBsInc,lBasl,lBsInc,iPrimi,iPrInc,jPrimj, &
            jPrInc,kPrimk,kPrInc,lPriml,lPrInc,ipMem1,ipMem2,Mem1,Mem2,MemPSO)
#ifdef _DEBUGPRINT_
write(u6,*) ' ************** Memory partioning **************'
write(u6,*) ' ipMem1=',ipMem1
write(u6,*) ' Mem1=',Mem1
write(u6,*) ' Mem2=',Mem2
write(u6,*) ' MemPSO=',MemPSO

write(u6,*) ' iBasi,iBsInc=',iBasi,iBsInc
write(u6,*) ' jBasj,jBsInc=',jBasj,jBsInc
write(u6,*) ' kBasi,kBsInc=',kBask,kBsInc
write(u6,*) ' lBasl,lBsInc=',lBasl,lBsInc
write(u6,*) ' iPrimi,iPrInc=',iPrimi,iPrInc
write(u6,*) ' jPrimj,jPrInc=',jPrimj,jPrInc
write(u6,*) ' kPrimk,kPrInc=',kPrimk,kPrInc
write(u6,*) ' lPriml,lPrInc=',lPriml,lPrInc
write(u6,*) ' ***********************************************'
#endif

return

end subroutine SOAO_g
