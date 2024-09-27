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
! Copyright (C) 1999, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine HrrMtrx(HMtrx,np,la,lb,A,B,Sph_a,CS_a,nSph_a,Sph_b,Cs_b,nSph_b)
!***********************************************************************
!                                                                      *
!     Object: to compute the matrix which corresponds to the transfer  *
!             equation.                                                *
!                                                                      *
!     Author: Roland Lindh                                             *
!             Dept of Chem. Phys.                                      *
!             Univ. of Lund, Sweden                                    *
!             February 1999                                            *
!***********************************************************************

use Index_Functions, only: C3_Ind3, nTri_Elem1, nTri3_Elem
use define_af, only: Binom, iCan, iTabMx
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: np, la, lb, nSph_a, nSph_b
real(kind=wp), intent(out) :: HMtrx(np,nSph_a,nSph_b)
real(kind=wp), intent(in) :: A(3), B(3), CS_a(nTri_Elem1(la),nSph_a), CS_b(nTri_Elem1(lb),nSph_b)
logical(kind=iwp), intent(in) :: Sph_a, Sph_b
integer(kind=iwp) :: i, ipa, ipb, ipe, iSph_a, iSph_b, ix, ixLow, iy, iyLow, iz, izLow, jOff, jx, jxLow, jy, jyLow, jz, jzLow, kx, &
                     ky, kz, lta, ltb, na, nab, nb
real(kind=wp) :: AB(3,0:iTabMx), ABx, ABy, ABz, C_A, C_B
logical(kind=iwp), external :: EQ
#ifdef _DEBUGPRINT_
real(kind=wp), external :: DDot_
#endif

lta = nTri_Elem1(la)
ltb = nTri_Elem1(lb)
#ifdef _DEBUGPRINT_
call RecPrt('A',' ',A,1,3)
call RecPrt('B',' ',B,1,3)
call RecPrt('CS_a',' ',CS_a,lta,nSph_a)
call RecPrt('CS_b',' ',CS_b,ltb,nSph_b)
write(u6,*) 'np=',np
#endif
HMtrx(:,:,:) = Zero

AB(:,0) = One
if (la >= lb) then
  AB(:,1) = A(:)-B(:)
else
  AB(:,1) = B(:)-A(:)
end if
do i=2,min(la,lb)
  AB(:,i) = AB(:,i-1)*AB(:,1)
end do

na = nTri3_Elem(la)
nb = nTri3_Elem(lb)
nab = nTri3_Elem(la+lb)

if (la >= lb) then

  if (Sph_a .and. Sph_b) then

    do iSph_a=1,nSph_a
      do ipa=1,lta
        C_a = CS_a(ipa,iSph_a)
        if (C_a == Zero) cycle
        ix = iCan(1,na+ipa)
        iy = iCan(2,na+ipa)
        iz = iCan(3,na+ipa)
        do iSph_b=1,nSph_b
          do ipb=1,ltb
            C_b = CS_b(ipb,iSph_b)
            if (C_b == Zero) cycle
            jx = iCan(1,nb+ipb)
            jy = iCan(2,nb+ipb)
            jz = iCan(3,nb+ipb)

            if (EQ(A,B)) then
              ixLow = ix+jx
              iyLow = iy+jy
              izLow = iz+jz
              jOff = nab
            else
              ixLow = ix
              iyLow = iy
              izLow = iz
              jOff = na
            end if
            do kx=ixLow,ix+jx
              do ky=iyLow,iy+jy
                do kz=izLow,iz+jz
                  ipe = C3_Ind3(kx,ky,kz)-jOff

                  ABx = AB(1,ix+jx-kx)*Binom(jx,kx-ix)
                  ABy = AB(2,iy+jy-ky)*Binom(jy,ky-iy)
                  ABz = AB(3,iz+jz-kz)*Binom(jz,kz-iz)
                  HMtrx(ipe,iSph_a,iSph_b) = HMtrx(ipe,iSph_a,iSph_b)+ABx*ABy*ABz*C_a*C_b
                end do
              end do
            end do

          end do
        end do
      end do
    end do

  else if (Sph_a) then

    do iSph_a=1,nSph_a
      do ipa=1,lta
        C_a = CS_a(ipa,iSph_a)
        if (C_a == Zero) cycle
        ix = iCan(1,na+ipa)
        iy = iCan(2,na+ipa)
        iz = iCan(3,na+ipa)
        do ipb=1,ltb
          jx = iCan(1,nb+ipb)
          jy = iCan(2,nb+ipb)
          jz = iCan(3,nb+ipb)

          if (EQ(A,B)) then
            ixLow = ix+jx
            iyLow = iy+jy
            izLow = iz+jz
            jOff = nab
          else
            ixLow = ix
            iyLow = iy
            izLow = iz
            jOff = na
          end if
          do kx=ixLow,ix+jx
            do ky=iyLow,iy+jy
              do kz=izLow,iz+jz
                ipe = C3_Ind3(kx,ky,kz)-jOff

                ABx = AB(1,ix+jx-kx)*Binom(jx,kx-ix)
                ABy = AB(2,iy+jy-ky)*Binom(jy,ky-iy)
                ABz = AB(3,iz+jz-kz)*Binom(jz,kz-iz)
                HMtrx(ipe,iSph_a,ipb) = HMtrx(ipe,iSph_a,ipb)+ABx*ABy*ABz*C_a
              end do
            end do
          end do

        end do
      end do
    end do

  else if (Sph_b) then

    do ipa=1,lta
      ix = iCan(1,na+ipa)
      iy = iCan(2,na+ipa)
      iz = iCan(3,na+ipa)
      do iSph_b=1,nSph_b
        do ipb=1,ltb
          C_b = CS_b(ipb,iSph_b)
          if (C_b == Zero) cycle
          jx = iCan(1,nb+ipb)
          jy = iCan(2,nb+ipb)
          jz = iCan(3,nb+ipb)

          if (EQ(A,B)) then
            ixLow = ix+jx
            iyLow = iy+jy
            izLow = iz+jz
            jOff = nab
          else
            ixLow = ix
            iyLow = iy
            izLow = iz
            jOff = na
          end if
          do kx=ixLow,ix+jx
            do ky=iyLow,iy+jy
              do kz=izLow,iz+jz
                ipe = C3_Ind3(kx,ky,kz)-jOff

                ABx = AB(1,ix+jx-kx)*Binom(jx,kx-ix)
                ABy = AB(2,iy+jy-ky)*Binom(jy,ky-iy)
                ABz = AB(3,iz+jz-kz)*Binom(jz,kz-iz)
                HMtrx(ipe,ipa,iSph_b) = HMtrx(ipe,ipa,iSph_b)+ABx*ABy*ABz*C_b
              end do
            end do
          end do

        end do
      end do
    end do

  else

    do ipa=1,lta
      ix = iCan(1,na+ipa)
      iy = iCan(2,na+ipa)
      iz = iCan(3,na+ipa)
      do ipb=1,ltb
        jx = iCan(1,nb+ipb)
        jy = iCan(2,nb+ipb)
        jz = iCan(3,nb+ipb)

        if (EQ(A,B)) then
          ixLow = ix+jx
          iyLow = iy+jy
          izLow = iz+jz
          jOff = nab
        else
          ixLow = ix
          iyLow = iy
          izLow = iz
          jOff = na
        end if
        do kx=ixLow,ix+jx
          do ky=iyLow,iy+jy
            do kz=izLow,iz+jz
              ipe = C3_Ind3(kx,ky,kz)-jOff

              ABx = AB(1,ix+jx-kx)*Binom(jx,kx-ix)
              ABy = AB(2,iy+jy-ky)*Binom(jy,ky-iy)
              ABz = AB(3,iz+jz-kz)*Binom(jz,kz-iz)
              HMtrx(ipe,ipa,ipb) = HMtrx(ipe,ipa,ipb)+ABx*ABy*ABz
            end do
          end do
        end do

      end do
    end do

  end if

else

  if (Sph_a .and. Sph_b) then

    do iSph_a=1,nSph_a
      do ipa=1,lta
        C_a = CS_a(ipa,iSph_a)
        if (C_a == Zero) cycle
        ix = iCan(1,na+ipa)
        iy = iCan(2,na+ipa)
        iz = iCan(3,na+ipa)
        do iSph_b=1,nSph_b
          do ipb=1,ltb
            C_b = CS_b(ipb,iSph_b)
            if (C_b == Zero) cycle
            jx = iCan(1,nb+ipb)
            jy = iCan(2,nb+ipb)
            jz = iCan(3,nb+ipb)

            if (EQ(A,B)) then
              jxLow = ix+jx
              jyLow = iy+jy
              jzLow = iz+jz
              jOff = nab
            else
              jxLow = jx
              jyLow = jy
              jzLow = jz
              jOff = nb
            end if
            do kx=jxLow,ix+jx
              do ky=jyLow,iy+jy
                do kz=jzLow,iz+jz
                  ipe = C3_Ind3(kx,ky,kz)-jOff

                  ABx = AB(1,ix+jx-kx)*Binom(ix,kx-jx)
                  ABy = AB(2,iy+jy-ky)*Binom(iy,ky-jy)
                  ABz = AB(3,iz+jz-kz)*Binom(iz,kz-jz)
                  HMtrx(ipe,iSph_a,iSph_b) = HMtrx(ipe,iSph_a,iSph_b)+ABx*ABy*ABz*C_a*C_b
                end do
              end do
            end do

          end do
        end do
      end do
    end do

  else if (Sph_a) then

    do iSph_a=1,nSph_a
      do ipa=1,lta
        C_a = CS_a(ipa,iSph_a)
        if (C_a == Zero) cycle
        ix = iCan(1,na+ipa)
        iy = iCan(2,na+ipa)
        iz = iCan(3,na+ipa)
        do ipb=1,ltb
          jx = iCan(1,nb+ipb)
          jy = iCan(2,nb+ipb)
          jz = iCan(3,nb+ipb)

          if (EQ(A,B)) then
            jxLow = ix+jx
            jyLow = iy+jy
            jzLow = iz+jz
            jOff = nab
          else
            jxLow = jx
            jyLow = jy
            jzLow = jz
            jOff = nb
          end if
          do kx=jxLow,ix+jx
            do ky=jyLow,iy+jy
              do kz=jzLow,iz+jz
                ipe = C3_Ind3(kx,ky,kz)-jOff

                ABx = AB(1,ix+jx-kx)*Binom(ix,kx-jx)
                ABy = AB(2,iy+jy-ky)*Binom(iy,ky-jy)
                ABz = AB(3,iz+jz-kz)*Binom(iz,kz-jz)
                HMtrx(ipe,iSph_a,ipb) = HMtrx(ipe,iSph_a,ipb)+ABx*ABy*ABz*C_a
              end do
            end do
          end do

        end do
      end do
    end do

  else if (Sph_b) then

    do ipa=1,lta
      ix = iCan(1,na+ipa)
      iy = iCan(2,na+ipa)
      iz = iCan(3,na+ipa)
      do iSph_b=1,nSph_b
        do ipb=1,ltb
          C_b = CS_b(ipb,iSph_b)
          if (C_b == Zero) cycle
          jx = iCan(1,nb+ipb)
          jy = iCan(2,nb+ipb)
          jz = iCan(3,nb+ipb)

          if (EQ(A,B)) then
            jxLow = ix+jx
            jyLow = iy+jy
            jzLow = iz+jz
            jOff = nab
          else
            jxLow = jx
            jyLow = jy
            jzLow = jz
            jOff = nb
          end if
          do kx=jxLow,ix+jx
            do ky=jyLow,iy+jy
              do kz=jzLow,iz+jz
                ipe = C3_Ind3(kx,ky,kz)-jOff

                ABx = AB(1,ix+jx-kx)*Binom(ix,kx-jx)
                ABy = AB(2,iy+jy-ky)*Binom(iy,ky-jy)
                ABz = AB(3,iz+jz-kz)*Binom(iz,kz-jz)
                HMtrx(ipe,ipa,iSph_b) = HMtrx(ipe,ipa,iSph_b)+ABx*ABy*ABz*C_b
              end do
            end do
          end do

        end do
      end do
    end do

  else

    do ipa=1,lta
      ix = iCan(1,na+ipa)
      iy = iCan(2,na+ipa)
      iz = iCan(3,na+ipa)
      do ipb=1,ltb
        jx = iCan(1,nb+ipb)
        jy = iCan(2,nb+ipb)
        jz = iCan(3,nb+ipb)

        if (EQ(A,B)) then
          jxLow = ix+jx
          jyLow = iy+jy
          jzLow = iz+jz
          jOff = nab
        else
          jxLow = jx
          jyLow = jy
          jzLow = jz
          jOff = nb
        end if
        do kx=jxLow,ix+jx
          do ky=jyLow,iy+jy
            do kz=jzLow,iz+jz
              ipe = C3_Ind3(kx,ky,kz)-jOff

              ABx = AB(1,ix+jx-kx)*Binom(ix,kx-jx)
              ABy = AB(2,iy+jy-ky)*Binom(iy,ky-jy)
              ABz = AB(3,iz+jz-kz)*Binom(iz,kz-jz)
              HMtrx(ipe,ipa,ipb) = HMtrx(ipe,ipa,ipb)+ABx*ABy*ABz
            end do
          end do
        end do

      end do
    end do

  end if

end if
#ifdef _DEBUGPRINT_
call RecPrt('HMat ( np x (nSph_a*nSph_b) )','(30F4.1)',HMtrx,np,nSph_a*nSph_b)
write(u6,*) DDot_(np*nSph_a*nSph_b,HMtrx,1,HMtrx,1)
#endif

return

end subroutine HrrMtrx
