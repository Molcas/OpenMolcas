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

use Index_Functions, only: C_Ind3, nTri_Elem1
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iAng, nCoor, nElem, nCmp, nTerm, nForm, nRad, mExp, nExp, nBas, mAO, ipx, ipy, ipz
real(kind=wp), intent(in) :: Coor(3,nCoor), RA(3), CffSph(nElem,nCmp), Thr_Rad, Alpha(nExp), CffCnt(mExp,nBas), px, py, pz
real(kind=wp), intent(out) :: xyz(nCoor,3,0:iAng+nRad-1), Radial(nCoor,nRad,nBas), AOValue(mAO,nCoor,nBas,nCmp)
logical(kind=iwp), intent(in) :: Transf
integer(kind=iwp), intent(out) :: Angular(nTerm,5,nForm)
integer(kind=iwp) :: iX, iY, iZ, nDrv, iExp, iCoor, iBas, iDrv, i, i_f, ip, jDrv, jX, jY, jZ, jF, kDrv, iCmp, kForm, mForm, mTerm, &
                     iTerm, iCoef, iRad, iForm
real(kind=wp) :: Cff, Coef, R2, ThrE, Tmp
#ifdef _DEBUGPRINT_
character(len=80) :: Label
#endif

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
Radial(:,:,:) = Zero
do iCoor=1,nCoor
  R2 = (Coor(1,iCoor)-RA(1))**2+(Coor(2,iCoor)-RA(2))**2+(Coor(3,iCoor)-RA(3))**2
  do iExp=1,nExp
    if (-Alpha(iExp)*R2 < Thre) exit
    Tmp = exp(-Alpha(iExp)*R2)
    Radial(iCoor,1,:) = Radial(iCoor,1,:)+CffCnt(iExp,:)*Tmp
    do iDrv=1,nDrv
      Tmp = -Two*Alpha(iExp)*Tmp
      Radial(iCoor,iDrv+1,:) = Radial(iCoor,iDrv+1,:)+CffCnt(iExp,:)*Tmp
    end do
  end do
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
  xyz(:,:,0) = One
  xyz(:,1,1) = px*(Coor(1,:)-RA(1))
  xyz(:,2,1) = py*(Coor(2,:)-RA(2))
  xyz(:,3,1) = pz*(Coor(3,:)-RA(3))
  do i=2,iAng+nRad-1
    xyz(:,:,i) = xyz(:,:,i-1)*xyz(:,:,1)
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
    ip = C_Ind3(ix,iy,iz)

    ! Initial values for 0th-order derivative

    Angular(1,1,1) = ix
    Angular(1,2,1) = iy
    Angular(1,3,1) = iz
    Angular(1,4,1) = 0
    Angular(1,5,1) = 1

    ! Calculate derivatives of the angular components

    ! jf: Formula to be derivated
    ! i_f: New formula by way of derivation

    i_f = 1
    !call PrntA(nterm,nform,Angular,i_f,0)
    do jDrv=0,nDrv-1
      do jx=jDrv,0,-1
        do jy=jDrv-jx,0,-1
          jz = jDrv-jx-jy
          jf = C_Ind3(jx,jy,jz)

          do kDrv=0,jDrv-1
            jf = jf+nTri_Elem1(kDrv)
          end do

#         ifdef _DEBUGPRINT_
          write(u6,*) ' jx,jy,jz,jf=',jx,jy,jz,jf
#         endif

          if ((jy == 0) .and. (jz == 0)) then
            i_f = i_f+1
            call dFdxyz(nTerm,nForm,Angular,jf,i_f,1,ipx,jDrv)
            !call PrntA(nterm,nform,Angular,i_f,jdrv+1)
            i_f = i_f+1
            call dFdxyz(nTerm,nForm,Angular,jf,i_f,2,ipy,jDrv)
            !call PrntA(nterm,nform,Angular,i_f,jdrv+1)
            i_f = i_f+1
            call dFdxyz(nTerm,nForm,Angular,jf,i_f,3,ipz,jDrv)
            !call PrntA(nterm,nform,Angular,i_f,jdrv+1)
          else if (jz == 0) then
            i_f = i_f+1
            call dFdxyz(nTerm,nForm,Angular,jf,i_f,2,ipy,jDrv)
            !call PrntA(nterm,nform,Angular,i_f,jdrv+1)
            i_f = i_f+1
            call dFdxyz(nTerm,nForm,Angular,jf,i_f,3,ipz,jDrv)
            !call PrntA(nterm,nform,Angular,i_f,jdrv+1)
          else
            i_f = i_f+1
            call dFdxyz(nTerm,nForm,Angular,jf,i_f,3,ipz,jDrv)
            !call PrntA(nterm,nform,Angular,i_f,jdrv+1)
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
            mForm = nTri_Elem1(iDrv)
            do iForm=kForm+1,kForm+mForm
              mTerm = 2**(iDrv)
              do iTerm=1,mTerm
                iCoef = Angular(iTerm,5,iForm)
                if (iCoef /= 0) then
                  Coef = real(iCoef,kind=wp)
                  iRad = Angular(iTerm,4,iForm)+1
                  do iBas=1,nBas
                    AOValue(iForm,:,iBas,iCmp) = AOValue(iForm,:,iBas,iCmp)+xyz(:,1,Angular(iTerm,1,iForm))* &
                                                 xyz(:,2,Angular(iTerm,2,iForm))*xyz(:,3,Angular(iTerm,3,iForm))* &
                                                 Coef*Cff*Radial(:,iRad,iBas)
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
        mForm = nTri_Elem1(iDrv)
        do iForm=kForm+1,kForm+mForm
          mTerm = 2**(iDrv)
          do iTerm=1,mTerm
            iCoef = Angular(iTerm,5,iForm)
            if (iCoef /= 0) then
              Coef = real(iCoef,kind=wp)
              iRad = Angular(iTerm,4,iForm)+1
              do iBas=1,nBas
                AOValue(iForm,:,iBas,ip) = AOValue(iForm,:,iBas,ip)+xyz(:,1,Angular(iTerm,1,iForm))* &
                                           xyz(:,2,Angular(iTerm,2,iForm))*xyz(:,3,Angular(iTerm,3,iForm))* &
                                           Coef*Radial(:,iRad,iBas)
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
