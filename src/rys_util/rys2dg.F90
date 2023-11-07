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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine Rys2Dg(xyz2D0,nT,nRys,la,lb,lc,ld,xyz2D1,IfGrad,IndGrd,Coora,Alpha,Beta,Gmma,Delta,nZeta,nEta,Scrtch,Temp,Indx,ExpX, &
                  ExpY,mZeta,mEta)
!***********************************************************************
!                                                                      *
! Object: to compute the gradients of the 2D-integrals.                *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use Constants, only: One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nT, nRys, la, lb, lc, ld, nZeta, nEta, mZeta, mEta
real(kind=wp), intent(in) :: xyz2D0(nRys*nT,0:la+1,0:lb+1,0:lc+1,0:ld+1,3), Coora(3,4), Alpha(nZeta), Beta(nZeta), Gmma(nEta), &
                             Delta(nEta)
real(kind=wp), intent(out) :: xyz2D1(nRys*nT,0:la,0:lb,0:lc,0:ld,3,3), Scrtch(nRys*nT), Temp(nT)
logical(kind=iwp), intent(inout) :: IfGrad(3,4)
integer(kind=iwp), intent(inout) :: IndGrd(3,4), Indx(3,4)
external :: ExpX, ExpY
integer(kind=iwp) :: i1, i2, ia, ib, ic, iCar, iCent, id, Ind1(3), Ind2(3), jCent, nVec, nx, ny, nz
real(kind=wp) :: Fact
logical(kind=iwp), external :: EQ

#ifdef _DEBUGPRINT_
call RecPrt(' In Rys2Dg: Alpha',' ',Alpha,1,nZeta)
call RecPrt(' In Rys2Dg: Beta ',' ',Beta,1,nZeta)
call RecPrt(' In Rys2Dg: Gamma',' ',Gmma,1,nEta)
call RecPrt(' In Rys2Dg: Delta',' ',Delta,1,nEta)
write(u6,*) ' IfGrad=',IfGrad
write(u6,*) ' IndGrd=',IndGrd
#endif
nx = 0
ny = 0
nz = 0
Indx(:,:) = 0

! Differentiate with respect to the first center

if (IfGrad(1,1) .or. IfGrad(2,1) .or. IfGrad(3,1)) then
  call ExpX(Temp,mZeta,mEta,Alpha,One)
  call Exp_2(Scrtch,nRys,nT,Temp,One)
  !call RecPrt('Expanded exponents (alpha)',' ',Scrtch,nT,nRys)
end if
nVec = 0
if (IfGrad(1,1)) then
  nx = nx+1
  nVec = nVec+1
  Ind1(nVec) = nx
  Ind2(nVec) = 1
  Indx(1,1) = nx
end if
if (IfGrad(2,1)) then
  ny = ny+1
  nVec = nVec+1
  Ind1(nVec) = ny
  Ind2(nVec) = 2
  Indx(2,1) = ny
end if
if (IfGrad(3,1)) then
  nz = nz+1
  nVec = nVec+1
  Ind1(nVec) = nz
  Ind2(nVec) = 3
  Indx(3,1) = nz
end if

if (nVec /= 0) then

  do id=0,ld
    do ic=0,lc
      do ib=0,lb
        if (nVec == 3) then
          xyz2D1(:,0,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,1,ib,ic,id,Ind2(1))
          xyz2D1(:,0,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,1,ib,ic,id,Ind2(2))
          xyz2D1(:,0,ib,ic,id,Ind2(3),Ind1(3)) = Two*Scrtch(:)*xyz2D0(:,1,ib,ic,id,Ind2(3))
          if (la >= 1) then
            Fact = -One
            do ia=1,la
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia+1,ib,ic,id,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia-1,ib,ic,id,Ind2(1))
              xyz2D1(:,ia,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia+1,ib,ic,id,Ind2(2))+ &
                                                      Fact*xyz2D0(:,ia-1,ib,ic,id,Ind2(2))
              xyz2D1(:,ia,ib,ic,id,Ind2(3),Ind1(3)) = Two*Scrtch(:)*xyz2D0(:,ia+1,ib,ic,id,Ind2(3))+ &
                                                      Fact*xyz2D0(:,ia-1,ib,ic,id,Ind2(3))
              Fact = Fact-One
            end do
          end if
        else if (nVec == 2) then
          xyz2D1(:,0,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,1,ib,ic,id,Ind2(1))
          xyz2D1(:,0,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,1,ib,ic,id,Ind2(2))
          if (la >= 1) then
            Fact = -One
            do ia=1,la
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia+1,ib,ic,id,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia-1,ib,ic,id,Ind2(1))
              xyz2D1(:,ia,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia+1,ib,ic,id,Ind2(2))+ &
                                                      Fact*xyz2D0(:,ia-1,ib,ic,id,Ind2(2))
              Fact = Fact-One
            end do
          end if
        else if (nVec == 1) then
          xyz2D1(:,0,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,1,ib,ic,id,Ind2(1))
          if (la >= 1) then
            Fact = -One
            do ia=1,la
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia+1,ib,ic,id,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia-1,ib,ic,id,Ind2(1))
              Fact = Fact-One
            end do
          end if
        end if
      end do
    end do
  end do

end if

! Differentiate with respect to the second center

if (IfGrad(1,2) .or. IfGrad(2,2) .or. IfGrad(3,2)) then
  call ExpX(Temp,mZeta,mEta,Beta,One)
  call Exp_2(Scrtch,nRys,nT,Temp,One)
  !call RecPrt('Expanded exponents (beta) ',' ',Scrtch,nT,nRys)
end if
nVec = 0
if (IfGrad(1,2)) then
  nx = nx+1
  nVec = nVec+1
  Ind1(nVec) = nx
  Ind2(nVec) = 1
  Indx(1,2) = nx
end if
if (IfGrad(2,2)) then
  ny = ny+1
  nVec = nVec+1
  Ind1(nVec) = ny
  Ind2(nVec) = 2
  Indx(2,2) = ny
end if
if (IfGrad(3,2)) then
  nz = nz+1
  nVec = nVec+1
  Ind1(nVec) = nz
  Ind2(nVec) = 3
  Indx(3,2) = nz
end if

if (nVec /= 0) then

  do id=0,ld
    do ic=0,lc
      do ia=0,la
        if (nVec == 3) then
          xyz2D1(:,ia,0,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,1,ic,id,Ind2(1))
          xyz2D1(:,ia,0,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,1,ic,id,Ind2(2))
          xyz2D1(:,ia,0,ic,id,Ind2(3),Ind1(3)) = Two*Scrtch(:)*xyz2D0(:,ia,1,ic,id,Ind2(3))
          if (lb >= 1) then
            Fact = -One
            do ib=1,lb
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib+1,ic,id,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia,ib-1,ic,id,Ind2(1))
              xyz2D1(:,ia,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib+1,ic,id,Ind2(2))+ &
                                                      Fact*xyz2D0(:,ia,ib-1,ic,id,Ind2(2))
              xyz2D1(:,ia,ib,ic,id,Ind2(3),Ind1(3)) = Two*Scrtch(:)*xyz2D0(:,ia,ib+1,ic,id,Ind2(3))+ &
                                                      Fact*xyz2D0(:,ia,ib-1,ic,id,Ind2(3))
              Fact = Fact-One
            end do
          end if
        else if (nVec == 2) then
          xyz2D1(:,ia,0,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,1,ic,id,Ind2(1))
          xyz2D1(:,ia,0,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,1,ic,id,Ind2(2))
          if (lb >= 1) then
            Fact = -One
            do ib=1,lb
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib+1,ic,id,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia,ib-1,ic,id,Ind2(1))
              xyz2D1(:,ia,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib+1,ic,id,Ind2(2))+ &
                                                      Fact*xyz2D0(:,ia,ib-1,ic,id,Ind2(2))
              Fact = Fact-One
            end do
          end if
        else if (nVec == 1) then
          xyz2D1(:,ia,0,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,1,ic,id,Ind2(1))
          if (lb >= 1) then
            Fact = -One
            do ib=1,lb
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib+1,ic,id,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia,ib-1,ic,id,Ind2(1))
              Fact = Fact-One
            end do
          end if
        end if
      end do
    end do
  end do

end if

! Differentiate with respect to the third center

if (IfGrad(1,3) .or. IfGrad(2,3) .or. IfGrad(3,3)) then
  call ExpY(Temp,mZeta,mEta,Gmma,One)
  call Exp_2(Scrtch,nRys,nT,Temp,One)
  !call RecPrt('Expanded exponents (Gamma)',' ',Scrtch,nT,nRys)
end if
nVec = 0
if (IfGrad(1,3)) then
  nx = nx+1
  nVec = nVec+1
  Ind1(nVec) = nx
  Ind2(nVec) = 1
  Indx(1,3) = nx
end if
if (IfGrad(2,3)) then
  ny = ny+1
  nVec = nVec+1
  Ind1(nVec) = ny
  Ind2(nVec) = 2
  Indx(2,3) = ny
end if
if (IfGrad(3,3)) then
  nz = nz+1
  nVec = nVec+1
  Ind1(nVec) = nz
  Ind2(nVec) = 3
  Indx(3,3) = nz
end if

if (nVec /= 0) then

  do id=0,ld
    do ib=0,lb
      do ia=0,la
        if (nVec == 3) then
          xyz2D1(:,ia,ib,0,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,1,id,Ind2(1))
          xyz2D1(:,ia,ib,0,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,1,id,Ind2(2))
          xyz2D1(:,ia,ib,0,id,Ind2(3),Ind1(3)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,1,id,Ind2(3))
          if (lc >= 1) then
            Fact = -One
            do ic=1,lc
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic+1,id,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic-1,id,Ind2(1))
              xyz2D1(:,ia,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic+1,id,Ind2(2))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic-1,id,Ind2(2))
              xyz2D1(:,ia,ib,ic,id,Ind2(3),Ind1(3)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic+1,id,Ind2(3))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic-1,id,Ind2(3))
              Fact = Fact-One
            end do
          end if
        else if (nVec == 2) then
          xyz2D1(:,ia,ib,0,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,1,id,Ind2(1))
          xyz2D1(:,ia,ib,0,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,1,id,Ind2(2))
          if (lc >= 1) then
            Fact = -One
            do ic=1,lc
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic+1,id,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic-1,id,Ind2(1))
              xyz2D1(:,ia,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic+1,id,Ind2(2))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic-1,id,Ind2(2))
              Fact = Fact-One
            end do
          end if
        else if (nVec == 1) then
          xyz2D1(:,ia,ib,0,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,1,id,Ind2(1))
          if (lc >= 1) then
            Fact = -One
            do ic=1,lc
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic+1,id,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic-1,id,Ind2(1))
              Fact = Fact-One
            end do
          end if
        end if
      end do
    end do
  end do

end if

! Differentiate with respect to the fourth center

if (IfGrad(1,4) .or. IfGrad(2,4) .or. IfGrad(3,4)) then
  call ExpY(Temp,mZeta,mEta,Delta,One)
  call Exp_2(Scrtch,nRys,nT,Temp,One)
  !call RecPrt('Expanded exponents (delta)',' ',Scrtch,nT,nRys)
end if
nVec = 0
if (IfGrad(1,4)) then
  nx = nx+1
  nVec = nVec+1
  Ind1(nVec) = nx
  Ind2(nVec) = 1
  Indx(1,4) = nx
end if
if (IfGrad(2,4)) then
  ny = ny+1
  nVec = nVec+1
  Ind1(nVec) = ny
  Ind2(nVec) = 2
  Indx(2,4) = ny
end if
if (IfGrad(3,4)) then
  nz = nz+1
  nVec = nVec+1
  Ind1(nVec) = nz
  Ind2(nVec) = 3
  Indx(3,4) = nz
end if

if (nVec /= 0) then

  do ic=0,lc
    do ib=0,lb
      do ia=0,la
        if (nVec == 3) then
          xyz2D1(:,ia,ib,ic,0,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,1,Ind2(1))
          xyz2D1(:,ia,ib,ic,0,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,1,Ind2(2))
          xyz2D1(:,ia,ib,ic,0,Ind2(3),Ind1(3)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,1,Ind2(3))
          if (ld >= 1) then
            Fact = -One
            do id=1,ld
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,id+1,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic,id-1,Ind2(1))
              xyz2D1(:,ia,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,id+1,Ind2(2))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic,id-1,Ind2(2))
              xyz2D1(:,ia,ib,ic,id,Ind2(3),Ind1(3)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,id+1,Ind2(3))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic,id-1,Ind2(3))
              Fact = Fact-One
            end do
          end if
        else if (nVec == 2) then
          xyz2D1(:,ia,ib,ic,0,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,1,Ind2(1))
          xyz2D1(:,ia,ib,ic,0,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,1,Ind2(2))
          if (ld >= 1) then
            Fact = -One
            do id=1,ld
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,id+1,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic,id-1,Ind2(1))
              xyz2D1(:,ia,ib,ic,id,Ind2(2),Ind1(2)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,id+1,Ind2(2))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic,id-1,Ind2(2))
              Fact = Fact-One
            end do
          end if
        else if (nVec == 1) then
          xyz2D1(:,ia,ib,ic,0,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,1,Ind2(1))
          if (ld >= 1) then
            Fact = -One
            do id=1,ld
              xyz2D1(:,ia,ib,ic,id,Ind2(1),Ind1(1)) = Two*Scrtch(:)*xyz2D0(:,ia,ib,ic,id+1,Ind2(1))+ &
                                                      Fact*xyz2D0(:,ia,ib,ic,id-1,Ind2(1))
              Fact = Fact-One
            end do
          end if
        end if
      end do
    end do
  end do

end if

! Sum over common centers

do iCent=1,3
  do jCent=iCent+1,4
    if (EQ(Coora(1,iCent),Coora(1,jCent))) then
      do iCar=1,3

        if (IfGrad(iCar,iCent) .and. IfGrad(iCar,jCent)) then
          ! Change flags so gradient will not be assembled and
          ! that there will be no contribution to the gradient.
          IfGrad(iCar,jCent) = .false.
          IndGrd(iCar,jCent) = 0
          i1 = Indx(iCar,iCent)
          i2 = Indx(iCar,jCent)
          xyz2D1(:,:,:,:,:,iCar,i1) = xyz2D1(:,:,:,:,:,iCar,i1)+xyz2D1(:,:,:,:,:,iCar,i2)
        end if

      end do
    end if
  end do
end do

!#ifdef _DEBUGPRINT_
!do iCn=1,4
!  do iCar=1,3
!    if (IfGrad(iCar,iCn)) then
!      ij = Indx(iCar,iCn)
!      do ia=0,la
!        do ib=0,lb
!          do ic=0,lc
!            do id=0,ld
!              write(Label,'(A,4(I2,'',''),A,'','',I2,A)') ' xyz2D1(',ia,ib,ic,id,ch(iCar),iCn,')'
!              if (iPrint >= 99) then
!                call RecPrt(Label,' ',xyz2d1(1,ia,ib,ic,id,iCar,ij),nT,nRys)
!              else
!                write(u6,'(A)') Label
!                write(u6,*) DDot_(nT*nRys,xyz2d1(1,ia,ib,ic,id,iCar,ij),1,xyz2d1(1,ia,ib,ic,id,iCar,ij),1)
!              end if
!            end do
!          end do
!        end do
!      end do
!    end if
!  end do
!end do
!#endif

return

end subroutine Rys2Dg
