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

use external_centers
use Phase_Info

implicit real*8(A-H,O-Z)
external TNAI, Fake, XCff2D, XRys2D
#include "itmax.fh"
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables
real*8 C(3), TC(3), Coori(3,4), CoorAC(3,2), ZFd((iTabMx+1)*(iTabMx+2)/2), ZRFd((iTabMx+1)*(iTabMx+2)/2)
logical EQ, NoLoop, NoSpecial
integer iAnga(4), iDCRT(0:7), iStb(0:7), jCoSet(8,8)
character ChOper(0:7)*3
data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
!nElem(ixyz) = 2*ixyz+1
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

iRout = 151
iPrint = nPrint(iRout)

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

! Loop over charges and dipole moments in the external field

nData = 3
do iOrdOp=0,nOrdOp

  iAnga(1) = la
  iAnga(2) = lb
  iAnga(3) = iOrdOp
  iAnga(4) = 0
  call dcopy_(3,A,1,Coori(1,1),1)
  call dcopy_(3,RB,1,Coori(1,2),1)
  mabMin = nabSz(max(la,lb)-1)+1
  mabMax = nabSz(la+lb)
  if (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1
  mcdMin = nabSz(iOrdOp-1)+1
  mcdMax = nabSz(iOrdOp)
  lab = (mabMax-mabMin+1)
  kab = nElem(la)*nElem(lb)
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
    call dcopy_(3,A,1,CoorAC(1,1),1)
  else
    call dcopy_(3,RB,1,CoorAC(1,1),1)
  end if

  ! Loop over centers of the external field.

  iDum = 0
  do iFd=1,nXF

    NoLoop = .true.
    do jElem=1,nElem(iOrdOp)
      ZFd(jElem) = XF(nData+jElem,iFd)
      ! Divide quadrupole diagonal by 2 due to different normalisation
      if ((iOrdOp == 2) .and. ((jElem == 1) .or. (jElem == 4) .or. (jElem == 6))) ZFd(jElem) = ZFd(jElem)*0.5d0
      NoLoop = NoLoop .and. (ZFd(jElem) == Zero)
    end do

    if (NoLoop) Go To 111
    ! Pick up the center coordinates
    C(1:3) = XF(1:3,iFd)

    if (iPrint >= 99) call RecPrt('C',' ',C,1,3)

    ! Generate stabilizer of C

    iChxyz = iChAtm(C)
    call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

    !-Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
    Fact = dble(nStabM)/dble(LmbdT)

    if (iPrint >= 99) then
      write(6,*) ' m      =',nStabM
      write(6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
      write(6,*) ' s      =',nStb
      write(6,'(9A)') '(S)=',(ChOper(iStb(ii)),ii=0,nStb-1)
      write(6,*) ' LambdaT=',LmbdT
      write(6,*) ' t      =',nDCRT
      write(6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
    end if

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)

      jElem = 0
      do ix=iOrdOp,0,-1
        if (mod(ix,2) == 0) then
          Factx = One
        else
          Factx = dble(iPhase(1,iDCRT(lDCRT)))
        end if
        do iy=iOrdOp-ix,0,-1
          if (mod(iy,2) == 0) then
            Facty = One
          else
            Facty = dble(iPhase(2,iDCRT(lDCRT)))
          end if
          iz = iOrdOp-ix-iy
          if (mod(iz,2) == 0) then
            Factz = One
          else
            Factz = dble(iPhase(3,iDCRT(lDCRT)))
          end if

          jElem = jElem+1
          ZRFd(jElem) = Factx*Facty*Factz*ZFd(jElem)
        end do
      end do

      call dcopy_(3,TC,1,CoorAC(1,2),1)
      call dcopy_(3,TC,1,Coori(1,3),1)
      call dcopy_(3,TC,1,Coori(1,4),1)

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

      do i=1,nElem(iOrdOp)
        if (ZRFd(i) /= Zero) call SymAdO(Array(ipI),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,-Fact*ZRFd(i))
        ipI = ipI+nZeta*nElem(la)*nElem(lb)
      end do

#     ifdef _DEBUGPRINT_
      write(6,*) (Fact*ZFd(i),i=1,nElem(iOrdOp))
      call RecPrt('Array(ip1)',' ',Array(ip1),nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2*nElem(iOrdOp))
      call RecPrt('Final',' ',final,nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2*nIC)
#     endif

    end do ! End loop over DCRs
111 continue
  end do   ! iFd

  nData = nData+nElem(iOrdOp)
end do     ! iOrdOp

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Alpha)
  call Unused_real_array(Beta)
  call Unused_integer(nHer)
  call Unused_real_array(CCoor)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine XFdInt
