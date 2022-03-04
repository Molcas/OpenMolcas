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
subroutine contone(L,oneoverR3,onecontr,Lmax,contcoeff,nprim,ncont,MxcontL,dummy,onecartx,onecartY,onecartZ,charge,oneonly)
!bs contracts one-electron integrals and multiplies with l,m-dependent
!bs factors for L-,L0,L+

implicit real*8(a-h,o-z)
dimension oneoverR3(*), onecontr(MxcontL,MxcontL,-Lmax:Lmax,3), contcoeff(nprim,ncont), dummy(ncont,ncont), &
          onecartx(MxcontL,MxcontL, (Lmax+Lmax+1)*(Lmax+1)), onecarty(MxcontL,MxcontL, (Lmax+Lmax+1)*(Lmax+1)), &
          onecartz(MxcontL,MxcontL, (Lmax+Lmax+1)*(Lmax+1))
logical oneonly
!Statement function
ipnt(I,J) = (max(i,j)*(max(i,j)-1))/2+min(i,j)

!bs first of all cleaning dummy and onecontr
do jrun=1,ncont
  do irun=1,ncont
    dummy(irun,jrun) = 0d0
  end do
end do
if (oneonly) then
  iprod = MxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)
  call dzero(onecartx,iprod)
  call dzero(onecarty,iprod)
  call dzero(onecartz,iprod)
end if
iprod = 3*(Lmax+lmax+1)*MxcontL*MxcontL
call dzero(onecontr,iprod)
!bs contract onto dummy
do icont2=1,ncont
  do icont1=1,ncont
    do iprim2=1,nprim
      do iprim1=1,nprim
        dummy(icont1,icont2) = dummy(icont1,icont2)+contcoeff(iprim1,icont1)*contcoeff(iprim2,icont2)*oneoverR3(ipnt(iprim1,iprim2))
      end do
    end do
  end do
end do
do icont2=1,ncont
  do icont1=1,ncont
    dummy(icont1,icont2) = dummy(icont1,icont2)*charge
  end do
end do
!bs start to add l,m dependent factors
do M=-L,L
  factormin = sqrt(dble(L*L-M*M+L+M))
  factor0 = dble(M)
  factorplus = sqrt(dble(L*L-M*M+L-M))
  do irun=1,ncont
    do jrun=1,ncont
      onecontr(irun,jrun,M,1) = dummy(jrun,irun)*factormin  ! L-minus
    end do
  end do
  do irun=1,ncont
    do jrun=1,ncont
      onecontr(irun,jrun,M,2) = dummy(jrun,irun)*factor0    ! L-0
    end do
  end do
  do irun=1,ncont
    do jrun=1,ncont
      onecontr(irun,jrun,M,3) = dummy(jrun,irun)*factorplus ! L-plus
    end do
  end do
end do
!bs make the final cartesian integrals
call cartoneX(L,Lmax,onecontr,ncont,MxcontL,onecartX(1,1,1))
call cartoneY(L,Lmax,onecontr,ncont,MxcontL,onecartY(1,1,1))
call cartoneZ(L,Lmax,onecontr,ncont,MxcontL,onecartZ(1,1,1))

return

end subroutine contone
