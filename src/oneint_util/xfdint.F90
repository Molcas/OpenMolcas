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
! Copyright (C) 1995, Roland Lindh                                     *
!***********************************************************************

subroutine XFdInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, April '95                               *
!***********************************************************************

use external_centers, only: nXF, XF
use Phase_Info, only: iPhase
use Index_Functions, only: nTri3_Elem1, nTri_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: i, iAnga(4), iChxyz, iDCRT(0:7), iDum, iFd, ii, iOrdOp, ip1, ip2, ip3, ipI, ipIn, iPrint, iRout, iStb(0:7), &
                     ix, iy, iz, jCoSet(8,8), jElem, kab, lab, labcd, lcd, lDCRT, LmbdT, mabMax, mabMin, mArr, mcdMax, mcdMin, &
                     nData, nDCRT, nFLOP, nMem, nOp, nStb, nT
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), Fact, Factx, Facty, Factz, TC(3)
logical(kind=iwp) :: NoLoop, NoSpecial
real(kind=wp), allocatable :: ZFd(:), ZRFd(:)
character(len=*), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: iChAtm, NrOpr
logical(kind=iwp), external :: EQ
external :: TNAI, Fake, XCff2D, XRys2D

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(nHer)
unused_var(CoorO)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 151
iPrint = nPrint(iRout)

rFinal(:,:,:,:) = Zero

call mma_allocate(ZFd,nTri_Elem1(nOrdOp),label='ZFd')
call mma_allocate(ZRFd,nTri_Elem1(nOrdOp),label='ZRFd')

! Loop over charges and dipole moments in the external field

nData = 3
do iOrdOp=0,nOrdOp

  iAnga(1) = la
  iAnga(2) = lb
  iAnga(3) = iOrdOp
  iAnga(4) = 0
  Coori(:,1) = A
  Coori(:,2) = RB
  mabMin = nTri3_Elem1(max(la,lb)-1)
  mabMax = nTri3_Elem1(la+lb)-1
  if (EQ(A,RB)) mabMin = nTri3_Elem1(la+lb-1)
  mcdMin = nTri3_Elem1(iOrdOp-1)
  mcdMax = nTri3_Elem1(iOrdOp)-1
  lab = (mabMax-mabMin+1)
  kab = nTri_Elem1(la)*nTri_Elem1(lb)
  lcd = (mcdMax-mcdMin+1)
  labcd = lab*lcd

  ! Compute FLOP's and size of work array which Hrr will use.

  call mHrr(la,lb,nFLOP,nMem)

  ! Distribute the work array

  ip2 = 1
  ip1 = ip2+nZeta*max(labcd,lcd*nMem)
  mArr = nArr-max(labcd,lcd*nMem)

  ! Find center to accumulate angular momentum on. (HRR)

  if (la >= lb) then
    CoorAC(:,1) = A
  else
    CoorAC(:,1) = RB
  end if

  ! Loop over centers of the external field.

  iDum = 0
  do iFd=1,nXF

    NoLoop = .true.
    do jElem=1,nTri_Elem1(iOrdOp)
      ZFd(jElem) = XF(nData+jElem,iFd)
      ! Divide quadrupole diagonal by 2 due to different normalisation
      if ((iOrdOp == 2) .and. ((jElem == 1) .or. (jElem == 4) .or. (jElem == 6))) ZFd(jElem) = ZFd(jElem)*Half
      NoLoop = NoLoop .and. (ZFd(jElem) == Zero)
    end do

    if (NoLoop) cycle
    ! Pick up the center coordinates
    C(1:3) = XF(1:3,iFd)

    if (iPrint >= 99) call RecPrt('C',' ',C,1,3)

    ! Generate stabilizer of C

    iChxyz = iChAtm(C)
    call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

    !-Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    if (iPrint >= 99) then
      write(u6,*) ' m      =',nStabM
      write(u6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
      write(u6,*) ' s      =',nStb
      write(u6,'(9A)') '(S)=',(ChOper(iStb(ii)),ii=0,nStb-1)
      write(u6,*) ' LambdaT=',LmbdT
      write(u6,*) ' t      =',nDCRT
      write(u6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
    end if

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)

      jElem = 0
      do ix=iOrdOp,0,-1
        if (mod(ix,2) == 0) then
          Factx = One
        else
          Factx = real(iPhase(1,iDCRT(lDCRT)),kind=wp)
        end if
        do iy=iOrdOp-ix,0,-1
          if (mod(iy,2) == 0) then
            Facty = One
          else
            Facty = real(iPhase(2,iDCRT(lDCRT)),kind=wp)
          end if
          iz = iOrdOp-ix-iy
          if (mod(iz,2) == 0) then
            Factz = One
          else
            Factz = real(iPhase(3,iDCRT(lDCRT)),kind=wp)
          end if

          jElem = jElem+1
          ZRFd(jElem) = Factx*Facty*Factz*ZFd(jElem)
        end do
      end do

      CoorAC(:,2) = TC
      Coori(:,3) = TC
      Coori(:,4) = TC

      ! Compute integrals with the Rys quadrature.

      nT = nZeta
      NoSpecial = .true.
      call Rys(iAnga,nT,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,rKappa,[One],Coori,Coori,CoorAC,mabmin,mabmax,mcdMin,mcdMax, &
               Array(ip1),mArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)

      ! The integrals are now ordered as ijkl,e,f

      ! a) Change the order to f,ijkl,e
      ! b) Unfold e to ab, f,ijkl,ab
      ! c) Change the order back to ijkl,ab,f

      ! a)

      call DGeTMO(Array(ip1),nZeta*lab,nZeta*lab,lcd,Array(ip2),lcd)

      ! b) Use the HRR to unfold e to ab

      call HRR(la,lb,A,RB,Array(ip2),lcd*nZeta,nMem,ipIn)
      ip3 = ip2-1+ipIn

      ! c)

      call DGeTMO(Array(ip3),lcd,lcd,nZeta*kab,Array(ip1),nZeta*kab)

      ! Accumulate contributions to the symmetry adapted operator

      nOp = NrOpr(iDCRT(lDCRT))
      ipI = ip1

      do i=1,nTri_Elem1(iOrdOp)
        if (ZRFd(i) /= Zero) call SymAdO(Array(ipI),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,-Fact*ZRFd(i))
        ipI = ipI+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)
      end do

#     ifdef _DEBUGPRINT_
      write(u6,*) (Fact*ZFd(i),i=1,nTri_Elem1(iOrdOp))
      call RecPrt('Array(ip1)',' ',Array(ip1),nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2*nTri_Elem1(iOrdOp))
      call RecPrt('rFinal',' ',rFinal,nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2*nIC)
#     endif

    end do ! End loop over DCRs
  end do   ! iFd

  nData = nData+nTri_Elem1(iOrdOp)
end do     ! iOrdOp

call mma_deallocate(ZFd)
call mma_deallocate(ZRFd)

return

end subroutine XFdInt
