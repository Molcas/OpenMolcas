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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RFNuc(CoOP,rNucMm,ir)
!***********************************************************************
!                                                                      *
! Object: to compute the multipole moments for the nuclei.             *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Basis_Info
use Center_Info
#ifdef _OBSOLETE_
use External_Centers, only: nOrd_XF, XF
#endif
use Phase_Info
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One

implicit none
integer ir
real*8 rNucMm((ir+1)*(ir+2)/2), CoOp(3)
real*8 A(3), RA(3)
integer iq, ix, iy, iz, ndc, iCnttp, iCnt, i, mdc
real*8 temp, ZA, CCoMx, CCoMy, CComZ
#ifdef _OBSOLETE_
real*8 rRMy(3)
integer iStb(0:7), jCoSet(0:7,0:7)
real*8 DAx, DAy, DAz, Qxx, Qxy, Qxz, Qyy, Qyz, Qzz
real*8 QRAxx, QRAyy, QRAzz, QRAxy, QRAxz, QRAyz
#endif
#ifdef _DEBUGPRINT_
integer ip
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt(' In RFNuc:CoOp',' ',CoOp,1,3)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the nuclear contribution to the multipole moments

! Contributions due to the charges of nuclear charges

iq = 0
do ix=ir,0,-1
  do iy=ir-ix,0,-1
    iq = iq+1
    iz = ir-ix-iy
    temp = Zero
    !write(6,*) ' ix,iy,iz=',ix,iy,iz

    ndc = 0
    do iCnttp=1,nCnttp
      ZA = dbsc(iCnttp)%Charge
      if (ZA == Zero) Go To 101
#     ifdef _DEBUGPRINT_
      write(6,*) ' Charge=',ZA
      call RecPrt(' Centers',' ',dbsc(iCnttp)%Coor,3,dbsc(iCnttp)%nCntr)
#     endif
      do iCnt=1,dbsc(iCnttp)%nCntr
        A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
        mdc = ndc+iCnt
        do i=0,nIrrep/dc(mdc)%nStab-1
          call OA(dc(mdc)%iCoSet(i,0),A,RA)
          !call RecPrt(' RA',' ',RA,1,3)
          !call RecPrt(' CoOp',' ',CoOp,1,3)
#         ifdef NAGFOR
          if (iCnt < -2) write(6,*) 'Nag problem'
#         endif

          if (ix == 0) then
            CCoMx = One
          else
            CCoMx = (RA(1)-CoOp(1))**ix
          end if
          if (iy == 0) then
            CCoMy = One
          else
            CCoMy = (RA(2)-CoOp(2))**iy
          end if
          if (iz == 0) then
            CCoMz = One
          else
            CCoMz = (RA(3)-CoOp(3))**iz
          end if
          !write(6,*) CCoMx, CCoMy, CCoMz, temp
          temp = temp+ZA*CCoMx*CCoMy*CCoMz
        end do
      end do
101   continue
      ndc = ndc+dbsc(iCnttp)%nCntr
    end do
    rNucMm(iq) = temp
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _OBSOLETE_
! The remainder of this subroutine is obsolete and
! is kept only for testing reasons. It is replaced by
! the subroutine XFMoment which is more general.
goto 99

if ((.not. allocated(XF)) .or. (nOrd_XF < 0)) Go To 99

! Contributions due to the charges and dipoles of the
! static external electric field.

!write(6,*) ' Adding contributions from esef!'

iq = 0
do ix=ir,0,-1
  do iy=ir-ix,0,-1
    iq = iq+1
    iz = ir-ix-iy
    temp = Zero
    !write(6,*) ' ix,iy,iz=',ix,iy,iz

    do iFd=1,nXF
      DAx = Zero
      DAy = Zero
      DAz = Zero
      Qxx = Zero
      Qxy = Zero
      Qxz = Zero
      Qyy = Zero
      Qyz = Zero
      Qzz = Zero
      if (nOrd_XF == 0) then
        ZA = XF(4,iFd)
      else if (nOrd_XF == 1) then
        ZA = XF(4,iFd)
        DAx = XF(5,iFd)
        DAy = XF(6,iFd)
        DAz = XF(7,iFd)
      else if (nOrd_XF == 2) then
        ZA = XF(4,iFd)
        DAx = XF(5,iFd)
        DAy = XF(6,iFd)
        DAz = XF(7,iFd)
        Qxx = XF(8,iFd)
        Qxy = XF(9,iFd)
        Qxz = XF(10,iFd)
        Qyy = XF(11,iFd)
        Qyz = XF(12,iFd)
        Qzz = XF(13,iFd)
      else
        call WarningMessage(2,'RFNuc: Option not implemented yet!')
        call Abend()
      end if
#     ifdef _DEBUGPRINT_
      write(6,*) ' Charge=',ZA
      write(6,*) ' ixyz=',ixyz
      call RecPrt(' Centers',' ',XF(1,iXF),3,1)
#     endif

      A(1:3) = XF(1:3,iXF)

      ! Generate Stabilazor of C

      iChxyz = iChAtm(A)
      iDum = 0
      call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

      !write(6,*) ' nStb=',nStb
      do i=0,nIrrep/nStb-1
        call OA(jCoSet(i,0),A,RA)
        rRMy(1) = DAx*dble(iPhase(1,jCoSet(i,0)))
        rRMy(2) = DAy*dble(iPhase(2,jCoSet(i,0)))
        rRMy(3) = DAz*dble(iPhase(3,jCoSet(i,0)))
        QRAxx = QAxx
        QRAyy = QAyy
        QRAzz = QAzz
        QRAxy = dble(iPhase(1,jCoSet(i,0))*iPhase(2,jCoSet(i,0)))*QAxy
        QRAxz = dble(iPhase(1,jCoSet(i,0))*iPhase(3,jCoSet(i,0)))*QAxz
        QRAyz = dble(iPhase(2,jCoSet(i,0))*iPhase(3,jCoSet(i,0)))*QAyz

        if (ix == 0) then
          CCoMx = One
        else
          CCoMx = (RA(1)-CoOp(1))**ix
        end if
        if (iy == 0) then
          CCoMy = One
        else
          CCoMy = (RA(2)-CoOp(2))**iy
        end if
        if (iz == 0) then
          CCoMz = One
        else
          CCoMz = (RA(3)-CoOp(3))**iz
        end if

        !write(6,*) CCoMx, CCoMy, CCoMz, temp

        ! The charge contribution

        temp = temp+ZA*CCoMx*CCoMy*CCoMz

        ! Dipole contributions

        if (ix >= 1) temp = temp+dble(ix)*rRmy(1)*CCoMy*CCoMz*(RA(1)-CoOp(1))**(ix-1)
        if (iy >= 1) temp = temp+dble(iy)*rRmy(2)*CCoMx*CCoMz*(RA(2)-CoOp(2))**(iy-1)
        if (iz >= 1) temp = temp+dble(iz)*rRmy(3)*CCoMx*CCoMy*(RA(3)-CoOp(3))**(iz-1)

      end do
    end do
    !write(6,*) ' Temp=',temp
    rNucMm(iq) = rNucMm(iq)+temp

  end do
end do

99 continue
#endif
#ifdef _DEBUGPRINT_
call RecPrt(' Nuclear Multipole Moments',' ',rNucMm,ip,1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine RFNuc
