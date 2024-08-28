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
!               2004, Takashi Tsuchiya                                 *
!***********************************************************************

subroutine AOEval(iAng,nCoor,Coor,xyz,RA,Transf,CffSph,nElem,nCmp,Angular,nTerm,nForm,Thr_Rad,nRad,mExp,nExp,Alpha,Radial,nBas, &
                  CffCnt,AOValue,mAO,px,py,pz,ipx,ipy,ipz)
!***********************************************************************
! Object: to compute the values of the AOs on a grid. The calculation  *
!         is sectioned in an angular part and a radial part. Due to the*
!         gradients we have 3 angular parts and 2 radial part.         *
!                                                                      *
!      Author:Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN. November '95                            *
!      Modified by: Takashi Tsuchiya, Dept. of Theoretical Chemistry,  *
!                   University of Lund, SWEDEN. February 2004          *
!***********************************************************************

use Constants, only: Zero, One, Two
use Definitions, only: wp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer nCoor, iAng, nRad, nElem, nCmp, nExp, nBas, mExp, nTerm, mAO, ipx, ipy, ipz, nForm
real*8 xyz(nCoor,3,0:iAng+nRad-1), Coor(3,nCoor), RA(3), CffSph(nElem,nCmp), Alpha(nExp), Radial(nCoor,nRad,nBas), &
       CffCnt(mExp,nBas), AOValue(mAO,nCoor,nBas,nCmp)
real*8 px, py, pz
integer Angular(nTerm,5,nForm)
logical Transf
#ifdef _DEBUGPRINT_
character(len=80) Label
#endif
integer iX, iY, iZ, nDrv, iExp, iCoor, iBas, iDrv, i, ip, if, jDrv, jX, jY, jZ, jF, kDrv, iCmp, kForm, mForm, mTerm, iTerm, iCoef, &
        iRad, Ind, iForm
real*8 ThrE, Thr_Rad, Exp_Min, R2, Tmp, XCff, Tmp2, Tmp3, Tmp4, Cff, Coef
!                                                                      *
!***********************************************************************
!                                                                      *
!     Index for derivatives in AOValue
!
!     mAO =   1:    F
!             2:   dFdx
!             3:   dFdy
!             4:   dFdz
!             5:  d2Fdxdx
!             6:  d2Fdxdy
!             7:  d2Fdxdz
!             8:  d2Fdydy
!             9:  d2Fdydz
!            10:  d2Fdzdz
!            11:  d3Fdxdxdx
!            12:  d3Fdxdxdy
!            13:  d3Fdxdxdz
!            14:  d3Fdxdydy
!            15:  d3Fdxdydz
!            16:  d3Fdxdzdz
!            17:  d3Fdydydy
!            18:  d3Fdydydz
!            19:  d3Fdydzdz
!            20:  d3Fdzdzdz
!            ... and so on;
!                in the same order for the higher order derivatives.
!                                                                      *
!***********************************************************************
!                                                                      *
! Statement function
Ind(iy,iz) = (iy+iz)*(iy+iz+1)/2+iz+1

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) '********** AOEval ***********'
write(u6,*) 'In AOEval'
call RecPrt('Coor',' ',Coor,3,nCoor)
call RecPrt('CffCnt',' ',CffCnt,mExp,nBas)
call RecPrt('CffSph',' ',CffSph,nElem,nCmp)
write(u6,*) 'RA=',RA
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize AOValue

AOValue(:,:,:,:) = Zero

!--- Set the order of derivation

nDrv = nRad-1
!write(u6,*) '----- nDrv = ', nDrv
!                                                                      *
!***********************************************************************
!                                                                      *
! Evaluate the radial part
!    Normal radial part and
!    premultiplied with (minus two times the exponent)**iDrv

#ifdef _DEBUGPRINT_
if ((nRad <= 0) .and. (nRad >= 5)) then
  write(u6,*) 'AOEval: illegal value of nRad!'
  call Abend()
end if
#endif

Thre = -99.0_wp
if (Thr_Rad > Zero) Thre = log(Thr_Rad)
Exp_Min = 1.0e10_wp
do iExp=1,nExp
  Exp_Min = min(Exp_Min,Alpha(iExp))
end do
Radial(:,:,:) = Zero
do iCoor=1,nCoor
  R2 = (Coor(1,iCoor)-RA(1))**2+(Coor(2,iCoor)-RA(2))**2+(Coor(3,iCoor)-RA(3))**2
  do iExp=1,nExp
    if (-Alpha(iExp)*R2 < Thre) Go To 9898
    Tmp = exp(-Alpha(iExp)*R2)
    if (nRad == 1) then
      do iBas=1,nBas
        XCff = CffCnt(iExp,iBas)
        Radial(iCoor,1,iBas) = Radial(iCoor,1,iBas)+XCff*Tmp
      end do
    else if (nRad == 2) then
      Tmp2 = -Two*Alpha(iExp)*Tmp
      do iBas=1,nBas
        XCff = CffCnt(iExp,iBas)
        Radial(iCoor,1,iBas) = Radial(iCoor,1,iBas)+XCff*Tmp
        Radial(iCoor,2,iBas) = Radial(iCoor,2,iBas)+XCff*Tmp2
      end do
    else if (nRad == 3) then
      Tmp2 = -Two*Alpha(iExp)*Tmp
      Tmp3 = -Two*Alpha(iExp)*Tmp2
      do iBas=1,nBas
        XCff = CffCnt(iExp,iBas)
        Radial(iCoor,1,iBas) = Radial(iCoor,1,iBas)+XCff*Tmp
        Radial(iCoor,2,iBas) = Radial(iCoor,2,iBas)+XCff*Tmp2
        Radial(iCoor,3,iBas) = Radial(iCoor,3,iBas)+XCff*Tmp3
      end do
    else if (nRad == 4) then
      Tmp2 = -Two*Alpha(iExp)*Tmp
      Tmp3 = -Two*Alpha(iExp)*Tmp2
      Tmp4 = -Two*Alpha(iExp)*Tmp3
      do iBas=1,nBas
        XCff = CffCnt(iExp,iBas)
        Radial(iCoor,1,iBas) = Radial(iCoor,1,iBas)+XCff*Tmp
        Radial(iCoor,2,iBas) = Radial(iCoor,2,iBas)+XCff*Tmp2
        Radial(iCoor,3,iBas) = Radial(iCoor,3,iBas)+XCff*Tmp3
        Radial(iCoor,4,iBas) = Radial(iCoor,4,iBas)+XCff*Tmp4
      end do
    else
      do iBas=1,nBas
        XCff = CffCnt(iExp,iBas)
        Radial(iCoor,1,iBas) = Radial(iCoor,1,iBas)+XCff*Tmp
      end do
      do iDrv=1,nDrv
        Tmp = -Two*Alpha(iExp)*Tmp
        do iBas=1,nBas
          XCff = CffCnt(iExp,iBas)
          Radial(iCoor,iDrv+1,iBas) = Radial(iCoor,iDrv+1,iBas)+XCff*Tmp
        end do
      end do
    end if
  end do
9898 continue
end do

#ifdef _DEBUGPRINT_
write(u6,*) mExp,nExp
write(Label,'(A)') 'Radial(nCoor*nRad,nBas)'
call RecPrt(Label,'(10G20.10)',Radial,nCoor*nRad,nBas)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Evaluate the angular part

if (iAng+nRad-1 >= 1) then
  do iCoor=1,nCoor
    xyz(iCoor,1,0) = One
    xyz(iCoor,2,0) = One
    xyz(iCoor,3,0) = One
    xyz(iCoor,1,1) = px*(Coor(1,iCoor)-RA(1))
    xyz(iCoor,2,1) = py*(Coor(2,iCoor)-RA(2))
    xyz(iCoor,3,1) = pz*(Coor(3,iCoor)-RA(3))
  end do
  do i=2,iAng+nRad-1
    do iCoor=1,nCoor
      xyz(iCoor,1,i) = xyz(iCoor,1,i-1)*xyz(iCoor,1,1)
      xyz(iCoor,2,i) = xyz(iCoor,2,i-1)*xyz(iCoor,2,1)
      xyz(iCoor,3,i) = xyz(iCoor,3,i-1)*xyz(iCoor,3,1)
    end do
  end do
else
  xyz(:,:,0) = One
end if
#ifdef _DEBUGPRINT_
do i=0,iAng+nRad-1
  write(Label,'(A,I2,A)') 'xyz(nCoor,nCar,',i,')'
  call RecPrt(Label,'(3ES13.6)',xyz(1,1,i),nCoor,3)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Calculate the angular components of the derivatives

Angular(:,:,:) = 0

do ix=iAng,0,-1

  do iy=iAng-ix,0,-1
    iz = iAng-ix-iy
    ip = Ind(iy,iz)

    ! Initial values for 0th-order derivative

    Angular(1,1,1) = ix
    Angular(1,2,1) = iy
    Angular(1,3,1) = iz
    Angular(1,4,1) = 0
    Angular(1,5,1) = 1

    ! Calculate derivatives of the angular components

    ! jf: Formula to be derivated
    ! if: New formula by way of derivation

    if = 1
    !call PrntA(nterm,nform,Angular,if,0)
    do jDrv=0,nDrv-1
      do jx=jDrv,0,-1
        do jy=jDrv-jx,0,-1
          jz = jDrv-jx-jy
          jf = Ind(jy,jz)

          do kDrv=0,jDrv-1
            jf = jf+(kDrv+1)*(kDrv+2)/2
          end do

#         ifdef _DEBUGPRINT_
          write(u6,*) ' jx,jy,jz,jf=',jx,jy,jz,jf
#         endif

          if ((jy == 0) .and. (jz == 0)) then
            if = if+1
            call dFdxyz(nTerm,nForm,Angular,jf,if,1,ipx,jDrv)
            !call PrntA(nterm,nform,Angular,if,jdrv+1)
            if = if+1
            call dFdxyz(nTerm,nForm,Angular,jf,if,2,ipy,jDrv)
            !call PrntA(nterm,nform,Angular,if,jdrv+1)
            if = if+1
            call dFdxyz(nTerm,nForm,Angular,jf,if,3,ipz,jDrv)
            !call PrntA(nterm,nform,Angular,if,jdrv+1)
          else if (jz == 0) then
            if = if+1
            call dFdxyz(nTerm,nForm,Angular,jf,if,2,ipy,jDrv)
            !call PrntA(nterm,nform,Angular,if,jdrv+1)
            if = if+1
            call dFdxyz(nTerm,nForm,Angular,jf,if,3,ipz,jDrv)
            !call PrntA(nterm,nform,Angular,if,jdrv+1)
          else
            if = if+1
            call dFdxyz(nTerm,nForm,Angular,jf,if,3,ipz,jDrv)
            !call PrntA(nterm,nform,Angular,if,jdrv+1)
          end if
        end do
      end do
    end do

    ! Distribute contributions to the real spherical harmonics
    !            and combine with angular part
    !       to yield the values and the gradients

    if (Transf) then
      do iCmp=1,nCmp
        Cff = CffSph(ip,iCmp)
        if (Cff /= Zero) then
          kForm = 0
          do iDrv=0,nDrv
            mForm = (iDrv+1)*(iDrv+2)/2
            do iForm=kForm+1,kForm+mForm
              mTerm = 2**(iDrv)
              do iTerm=1,mTerm
                iCoef = Angular(iTerm,5,iForm)
                if (iCoef /= 0) then
                  Coef = real(iCoef,kind=wp)
                  iRad = Angular(iTerm,4,iForm)+1
                  do iBas=1,nBas
                    do iCoor=1,nCoor
                      AOValue(iForm,iCoor,iBas,iCmp) = AOValue(iForm,iCoor,iBas,iCmp)+xyz(iCoor,1,Angular(iTerm,1,iForm))* &
                                                       xyz(iCoor,2,Angular(iTerm,2,iForm))*xyz(iCoor,3,Angular(iTerm,3,iForm))* &
                                                       Coef*Cff*Radial(iCoor,iRad,iBas)
                    end do
                  end do
                end if
              end do
            end do
            kForm = kForm+mForm
          end do
        end if
      end do
    else
      kForm = 0
      do iDrv=0,nDrv
        mForm = (iDrv+1)*(iDrv+2)/2
        do iForm=kForm+1,kForm+mForm
          mTerm = 2**(iDrv)
          do iTerm=1,mTerm
            iCoef = Angular(iTerm,5,iForm)
            if (iCoef /= 0) then
              Coef = real(iCoef,kind=wp)
              iRad = Angular(iTerm,4,iForm)+1
              do iBas=1,nBas
                do iCoor=1,nCoor
                  AOValue(iForm,iCoor,iBas,ip) = AOValue(iForm,iCoor,iBas,ip)+xyz(iCoor,1,Angular(iTerm,1,iForm))* &
                                                 xyz(iCoor,2,Angular(iTerm,2,iForm))*xyz(iCoor,3,Angular(iTerm,3,iForm))* &
                                                 Coef*Radial(iCoor,iRad,iBas)
                end do
              end do
            end if
          end do
        end do
        kForm = kForm+mForm
      end do
    end if

  end do

end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) 'mAO,nCoor,nBas,nCmp',mAO,nCoor,nBas,nCmp
call RecPrt('AOValue','(10G20.10)',AOValue,mAO,nCoor*nBas*nCmp)
#endif

end subroutine AOEval
