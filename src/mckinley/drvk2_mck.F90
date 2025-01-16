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

Subroutine Drvk2_mck()
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
use k2_arrays, only: DoGrad_, DoHess_
use k2_structure, only: k2data
use iSD_data, only: iSD, nSD
use Basis_Info, only: dbsc, Shells
use Symmetry_Info, only: iOper, nIrrep
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iAng, iBas, iCmpV(4), iCnt, iCnttp, iDCRR(0:7), ijCmp, ik2, ipM001, ipM002, &
                     ipM003, ipM004, iPrim, iPrInc, iS, iSD4(0:nSD,4), iShll, jAng, jBas, jCnt, jCnttp, &
                     jPrim, jPrInc, jS, jShll, kPrInc, lPrInc, M001, M002, M003, M004, M00d, MaxMem, MemPrm, &
                     nDCRR, nHrrab, nMemab, nSkal, nSO
real(kind=wp) :: Coor(3,2)
real(kind=wp), allocatable :: Wrk(:)

DoGrad_ = .false.
DoHess_ = .true.
call Nr_Shells(nSkal)
call Allok2()
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_MaxDBLE(MaxMem)
call mma_allocate(Wrk,(9*MaxMem)/10,Label='Wrk')
ipM001 = 1
!                                                                      *
!***********************************************************************
!                                                                      *
do iS=1,nSkal
  iSD4(:,1) = iSD(:,iS)
  iSD4(:,3) = iSD(:,iS)

  iShll = iSD4(0,1)
  iAng = iSD4(1,1)
  iBas = iSD4(3,1)
  iPrim = iSD4(5,1)
  iCnttp = iSD4(13,1)
  iCnt = iSD4(14,1)
  Coor(1:3,1) = dbsc(iCnttp)%Coor(1:3,iCnt)

  iCmpV(1) = nTri_Elem1(iAng)
  iSD4(2,1) = iCmpV(1)
  iSD4(2,3) = iCmpV(1)

  do jS=1,iS
    iSD4(:,2) = iSD(:,jS)
    iSD4(:,4) = iSD(:,jS)

    jShll = iSD4(0,2)
    jAng = iSD4(1,2)
    jBas = iSD4(3,2)
    jPrim = iSD4(5,2)
    jCnttp = iSD4(13,2)
    jCnt = iSD4(14,2)
    Coor(1:3,2) = dbsc(jCnttp)%Coor(1:3,jCnt)

    ik2 = iTri(iS,jS)
    iCmpV(2) = nTri_Elem1(jAng)
    iSD4(2,2) = iCmpV(2)
    iSD4(2,4) = iCmpV(2)

    ! Compute FLOP's for the transfer equation.

    call mHrr(iAng,jAng,nHrrab,nMemab)
    ijCmp = nTri_Elem1(iAng)*nTri_Elem1(jAng)

    iSD4(5,1) = iPrim
    iSD4(5,2) = jPrim

    iSD4(5,3) = 1
    iSD4(5,4) = 1

    iSD4(3,1) = iPrim
    iSD4(3,2) = jPrim
    iSD4(3,3) = 1
    iSD4(3,4) = 1

    iCmpV(3:4) = iCmpV(1:2)

    nSO = 1

    ! Compute memory request for the primitives, i.e. how much memory
    ! is needed up to the transfer equation.

    call MemRys(iSD4(1,:),MemPrm)

    ! Decide on the partioning of the shells based on
    ! the available memory and the requested memory.

    call PSOAO0_h(nSO,nMemab,nMemab,MemPrm,MaxMem,iPrInc,jPrInc,kPrInc,lPrInc,ipM001,ipM002,ipM003,ipM004,M001,M002,M003,M004, &
                  M00d,nSD,iSD4)
    if ((iSD4(3,1) /= iSD4(4,1)) .or. (iSD4(3,2) /= iSD4(4,2))) then
      write(u6,*) 'Drvk2: (iBasi /= iBsInc) .or. (jBasj /= jBsInc)'
      write(u6,*) 'iBasi,iBsInc=',iSD4(3,1),iSD4(4,1)
      write(u6,*) 'jBasj,jBsInc=',iSD4(3,2),iSD4(4,2)
      call Abend()
    end if

    ! Find the Double Coset Representatives for center A and B.

    iDCRR(0:nIrrep-1) = iOper(0:nIrrep-1)
    nDCRR = nIrrep

    ! Compute all pair entities (zeta, kappa, Px, Py, Pz, ZInv, alpha,
    ! beta, [nm|nm] and derivative entity, a total of ten different
    ! entities) for all possible unique pairs of centers generated
    ! for the symmetry unique centers A and B.

    call k2Loop_mck(Coor,iSD4(1,:),iDCRR,nDCRR,k2Data(:,ik2),ijCmp,Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim, &
                    Shells(iShll)%pCff,iBas,Shells(jShll)%pCff,jBas,nMemab,Wrk(ipM002),M002,Wrk(ipM003),M003)

  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Wrk)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Drvk2_mck
