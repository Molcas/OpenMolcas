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
! Copyright (C) 1990-1992, Roland Lindh                                *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Drvk2_mck(New_Fock)
!***********************************************************************
!                                                                      *
!  Object: to precompute all pair entites as zeta, kappa, P.           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             June '91, modified for k2 loop.                          *
!             January '92, modified to gradient calculations.          *
!             April '92, modified to use the Cauchy-Schwarz inequality *
!              to estimate the integral derivatives.                   *
!              Modified 1995 for 2nd derivatives by AB                 *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem1
use k2_arrays, only: DoGrad_, DoHess_, nDeDe
use k2_structure, only: Indk2, k2data
use iSD_data, only: iSD
use Basis_Info, only: dbsc, Shells
use Symmetry_Info, only: iOper, nIrrep
use Sizes_of_Seward, only: S
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: New_Fock
integer(kind=iwp) :: iAng, iAngV(4), iAO, iBas, iBasi, iBsInc, iCmp, iCmpV(4), iCnt, iCnttp, iDCRR(0:7), iDeSiz, ijCmp, ijShll, &
                     ik2, iShllV(2), ipM001, ipM002, ipM003, ipM004, iPrim, iPrimi, iPrInc, iS, iShell, iShll, iSmLbl, jAng, jAO, &
                     jBas, jBasj, jBsInc, jCmp, jCnt, jCnttp, jPrim, jPrimj, jPrInc, jS, jShell, jShll, kBask, kBsInc, kPrimk, &
                     kPrInc, lBasl, lBsInc, lPriml, lPrInc, M001, M002, M003, M004, M00d, MaxMem, MemPrm, MemTmp, mk2, nBasi, &
                     nBasj, nDCRR, nHrrab, nMemab, nSkal, nSO
real(kind=wp) :: Coor(3,2), TCpu1, TCpu2, TWall1, TWall2
real(kind=wp), allocatable :: Con(:), Wrk(:)
integer(kind=iwp), parameter :: nHm = 0
integer(kind=iwp), external :: MemSO1

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)

DoGrad_ = .false.
DoHess_ = .true.
call Nr_Shells(nSkal)
call Allok2()
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Con,S%m2Max,Label='Con')
!                                                                      *
!***********************************************************************
!                                                                      *
MemTmp = 0
do iAng=0,S%iAngMx
  MemTmp = max(MemTmp,(S%MaxPrm(iAng)*nTri_Elem1(iAng))**2)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_MaxDBLE(MaxMem)
call mma_allocate(Wrk,(9*MaxMem)/10,Label='Wrk')
ipM001 = 1
!                                                                      *
!***********************************************************************
!                                                                      *
ndede = 0
mk2 = 0
do iS=1,nSkal
  iShll = iSD(0,iS)
  iAng = iSD(1,iS)
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iPrim = iSD(5,iS)
  iAO = iSD(7,iS)
  iShell = iSD(11,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  Coor(1:3,1) = dbsc(iCnttp)%Coor(1:3,iCnt)

  iAngV(1) = iAng
  iShllV(1) = iShll
  iCmpV(1) = nTri_Elem1(iAng)

  do jS=1,iS
    jShll = iSD(0,jS)
    jAng = iSD(1,jS)
    jCmp = iSD(2,jS)
    jBas = iSD(3,jS)
    jPrim = iSD(5,jS)
    jAO = iSD(7,jS)
    jShell = iSD(11,jS)
    jCnttp = iSD(13,jS)
    jCnt = iSD(14,jS)
    Coor(1:3,2) = dbsc(jCnttp)%Coor(1:3,jCnt)

    ik2 = iTri(iS,jS)
    iAngV(2) = jAng
    iShllV(2) = jShll
    iCmpV(2) = nTri_Elem1(jAng)

    ! Compute FLOP's for the transfer equation.

    call mHrr(iAng,jAng,nHrrab,nMemab)
    ijCmp = nTri_Elem1(iAng)*nTri_Elem1(jAng)

    iPrimi = iPrim
    jPrimj = jPrim
    nBasi = Shells(iShllV(1))%nBasis
    nBasj = Shells(iShllV(2))%nBasis

    kPrimk = 1
    lPriml = 1
    iBasi = iPrimi
    jBasj = jPrimj
    kBask = 1
    lBasl = 1

    call ConMax(Con,iPrimi,jPrimj,Shells(iShll)%pCff,nBasi,Shells(jShll)%pCff,nBasj)

    iAngV(3:4) = iAngV(1:2)
    iCmpV(3:4) = iCmpV(1:2)

    ijShll = iTri(iShell,jShell)

    nSO = 1

    ! Compute memory request for the primitives, i.e. how much memory
    ! is needed up to the transfer equation.

    call MemRys(iAngV,MemPrm)

    ! Decide on the partioning of the shells based on
    ! the available memory and the requested memory.

    call PSOAO0_h(nSO,nMemab,nMemab,MemPrm,MaxMem,iAngV,iCmpV,iBasi,iBsInc,jBasj,jBsInc,kBask,kBsInc,lBasl,lBsInc,iPrimi,iPrInc, &
                  jPrimj,jPrInc,kPrimk,kPrInc,lPriml,lPrInc,ipM001,ipM002,ipM003,ipM004,M001,M002,M003,M004,M00d)
    if ((iBasi /= iBsInc) .or. (jBasj /= jBsInc)) then
      write(u6,*) 'Drvk2: (iBasi /= iBsInc) .or. (jBasj /= jBsInc)'
      write(u6,*) 'iBasi,iBsInc=',iBasi,iBsInc
      write(u6,*) 'jBasj,jBsInc=',jBasj,jBsInc
      call Abend()
    end if

    ! Find the Double Coset Representatives for center A and B.

    iDCRR(0:nIrrep-1) = iOper(0:nIrrep-1)
    nDCRR = nIrrep

    ! Compute all pair entities (zeta, kappa, Px, Py, Pz, ZInv, alpha,
    ! beta, [nm|nm] and derivative entity, a total of ten different
    ! entities) for all possible unique pairs of centers generated
    ! for the symmetry unique centers A and B.

    call k2Loop_mck(Coor,iAngV,iCmpV,iDCRR,nDCRR,k2Data(:,ik2), &
                    ijCmp,Shells(iShllV(1))%Exp,iPrimi,Shells(iShllV(2))%Exp,jPrimj, &
                    Shells(iShllV(1))%pCff,iBas,Shells(iShllV(2))%pCff,jBas,nMemab,Wrk(ipM002),M002,Wrk(ipM003),M003)

    Indk2(2,ijShll) = nDCRR
    Indk2(3,ijShll) = ik2
    mk2 = mk2+nDCRR

    if (New_Fock) then
      iDeSiz = 1+iPrim*jPrim+iCmp*jCmp
    else
      iDeSiz = 1+iPrim*jPrim+(iBas*jBas+1)*iCmp*jCmp
    end if
    iSmLbl = 1
    nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
    if (nSO > 0) nDeDe = nDeDe+iDeSiz*nDCRR

  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Wrk)
call mma_deallocate(Con)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,'(20X,A)') ' *** The k2 entities have been precomputed ***'
write(u6,'(I7,A)') mk2,' blocks of k2 data were computed.'
write(u6,'(A)') ' The prescreening is based on the integral estimates.'
#endif

call CWTime(TCpu2,TWall2)

return

end subroutine Drvk2_mck
