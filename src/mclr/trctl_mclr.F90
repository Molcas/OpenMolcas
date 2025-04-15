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

subroutine TRCTL_MCLR()
! Two-electron integral transformation program: control section
!
! Purpose: Set up of memory locations (decide if out of core is
!          needed. Loop over symmetry blocks of AO integrals.
!          The transformation routine TRAMO is called for each
!          symmetry block of integrals.

use Symmetry_Info, only: Mul
use MCLR_Data, only: CMO, FnHlf2, FnHlf3, FnTri1, FnTri2, FnTri3, FnTri4, FnTri5, ipCM, LuHlf2, LuHlf3, LuTri1, LuTri2, LuTri3, &
                     LuTri4, LuTri5
use input_mclr, only: nAsh, nBas, nFro, nIsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iAd13, iAD14, iAd23, iAd24, iAd34, iBatch, IntBuf, ipB, ipi, iSP, iSQ, iSR, iSS, lW1, lW2, lW3, lW4, lW5, &
                     MemX, nAP, nAQ, nAR, nAS, nBP, nBPQRS, nBQ, nBR, nBS, nDP, nDQ, nDR, nDS, nORBP, nSPQ, nSPQR, nSPQRS, nW1, &
                     nW2, nW3, nW4, nW5
integer(kind=iwp), allocatable :: Hlf1(:,:)
real(kind=wp), allocatable :: Buffer(:)
integer(kind=iwp), parameter :: LIOTAB = 512*512

!                                                                      *
!***********************************************************************
!                                                                      *

call DANAME_wa(LUTRI1,FNTRI1)
call DANAME_wa(LUTRI2,FNTRI2)
call DANAME_wa(LUTRI3,FNTRI3)
call DANAME_wa(LUTRI4,FNTRI4)
call DANAME_wa(LUTRI5,FNTRI5)
IAD14 = 0
IAD13 = 0
IAD23 = 0
IAD24 = 0
IAD34 = 0
call mma_allocate(Hlf1,LIOTAB,4,Label='Hlf1')

! Precompute start points

! Loop over quadruples of symmetries (nsp,nsq,nsr,nss)
! Note that the integrals on LUTWOAO have to be sorted in the
! same order as the loop structure below.

call mma_maxDBLE(MEMX)
MEMX = max(MEMX-MEMX/10,0)
call mma_allocate(Buffer,MEMX,Label='Buffer')
ipB = 1
IBATCH = 0
do iSP=1,NSYM
  NBP = NBAS(iSP)
  NAP = NASH(iSP)+nish(isp)
  nDP = nFro(iSP)!+nDel(iSP)
  do iSQ=1,iSP
    NBQ = NBAS(iSQ)
    NAQ = NASH(iSQ)+nish(isq)
    nDQ = nFro(iSQ)!+nDel(iSQ)
    NSPQ = Mul(iSP,iSQ)
    do iSR=1,NSYM
      NBR = NBAS(iSR)
      NAR = NASH(iSR)+nish(isr)
      nDR = nFro(iSR)!+nDel(iSR)
      NSPQR = Mul(NSPQ,iSR)
      do iSS=1,iSR
        NBS = NBAS(iSS)
        NAS = NASH(iSS)+nIsh(iSS)
        nDS = nFro(iSS)!+nDel(iSS)
        NSPQRS = Mul(NSPQR,iSS)

        ! Check the loop conditions and skip transformation step if possible

        NORBP = NAP*NAQ+NAR*NAS+NAP*NAR+NAP*NAS+NAQ*NAR+NAQ*NAS
        NBPQRS = NBP*NBQ*NBR*NBS
        if (NSPQRS /= 1) cycle
        IBATCH = IBATCH+1
        if (NORBP == 0) cycle
        if (NBPQRS == 0) cycle
        if (NAR+NAS == 0) cycle

        ! Set up dynamic memory

        INTBUF = 256*256
        ipi = ipB
        LW1 = ipi
        NW1 = max(NBP*NBQ,NBR*NBS)
        ipi = ipi+nw1
        LW2 = ipi
        NW2 = max(nBR*nAS,nAQ*nBP)
        ipi = ipi+nw2
        lw3 = ipi
        NW3 = max(nAR*nBS,nBQ*nAP)
        ipi = ipi+nw3
        lw4 = ipi
        NW4 = NAR*NAS
        ipi = ipi+nw4
        lw5 = ipi
        NW5 = MEMX-nw1-nw2-nw3-nw4

        ! transform the symmetry block (ISP,ISQ|ISR,ISS)

        call TRAMO_MCLR(INTBUF,Buffer(LW1:LW1+NW1-1),NW1,Buffer(LW2:LW2+NW2-1),NW2,Buffer(LW3:LW3+NW3-1),NW3, &
                        Buffer(LW4:LW4+NW4-1),NW4,Buffer(LW5:LW5+NW5-1),NW5,nBP,nBQ,nBR,nBS,iSP,iSQ,iSR,iSS,nAP,nAQ,nAR,nAS, &
                        CMO(ipCM(iSP)+nBP*nDP),CMO(ipCM(iSQ)+nBQ*nDQ),CMO(ipCM(iSR)+nBR*nDR),CMO(ipCM(iSS)+nBS*nDS),iAD13,iAD14, &
                        iAD23,iAD24,iAD34,Hlf1(:,1),Hlf1(:,2),Hlf1(:,3),Hlf1(:,4),LIOTAB)

        ! End of loop over quadruples of symmetries
      end do
    end do
  end do
end do

call mma_deallocate(Buffer)
call mma_deallocate(Hlf1)
call DAName_wa(LUHLF2,FNHLF2)
call DAName_wa(LUHLF3,FNHLF3)
call DaEras(LUHLF2)
call DaEras(LUHLF3)

end subroutine TRCTL_MCLR
