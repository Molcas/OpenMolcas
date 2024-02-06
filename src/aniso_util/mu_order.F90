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

subroutine mu_order(dim,MS,MM,gtens,order,HCF2,AMM,AMS,Z,iprint)
! This Subroutine receives the moment matrix dipso(3,dim,dim) and Returns the matrix re-builted using only the 1-st order operators.

use Constants, only: Zero, One, Half, cZero, Onei
use Definitions, only: u6

implicit none
integer, intent(in) :: dim, order, iprint
real(kind=8), intent(out) :: gtens(3)
! initial magnetic moment
complex(kind=8), intent(in) :: MM(3,dim,dim)
! initial spin moment
complex(kind=8), intent(in) :: MS(3,dim,dim)
! transformed magnetic moment
complex(kind=8), intent(out) :: AMM(3,dim,dim)
! transformed spin moment
complex(kind=8), intent(out) :: AMS(3,dim,dim)
complex(kind=8), intent(out) :: Z(dim,dim)
complex(kind=8), intent(out) :: HCF2(dim,3,dim,dim)
! local variables:
integer :: i, j, l, i1, i2, m, n
real(kind=8) :: maxes(3,3), m_fact
complex(kind=8) :: DIP_O(dim,dim), DIP_W(dim,dim), B(3,dim,-dim:dim), BNMC(3,dim,0:dim), BNMS(3,dim,0:dim), SP_MOW, SP_DIPO(3), &
                   SP_DIPW(3), O1, O2, trace
external :: trace
!------------------------------------------------------------

Z(:,:) = cZero
! get the local pseudospin:
call pseudospin(MM,dim,Z,3,1,iprint)
! re-write MM and MS to the new pseudospin basis:
AMS(:,:,:) = cZero
AMM(:,:,:) = cZero
do l=1,3
  do i=1,dim
    do j=1,dim
      do i1=1,dim
        do i2=1,dim
          AMM(l,i,j) = AMM(l,i,j)+MM(l,i1,i2)*conjg(Z(i1,i))*Z(i2,j)
          AMS(l,i,j) = AMS(l,i,j)+MS(l,i1,i2)*conjg(Z(i1,i))*Z(i2,j)
        end do
      end do
    end do
  end do
end do
if (iprint > 2) then
  call prMom('MU_ORDER:   AMM(l,i,j):',AMM,dim)
  call prMom('MU_ORDER:   AMS(l,i,j):',AMS,dim)
end if
! project the moment on ITO:
! obtain the b3m and c3m coefficients:
B(1:3,1:dim,-dim:dim) = cZero
do N=1,dim-1
  do M=0,N
    DIP_O(:,:) = cZero
    DIP_W(:,:) = cZero
    call Stewens_matrixel(N,M,dim,DIP_O,DIP_W,iprint)
    if (iprint > 5) then
      write(u6,*)
      write(u6,'( 5x,a,i2,a,i3)') 'DIP_O, n = ',N,', m =',m
      write(u6,*)
      do i=1,dim
        write(u6,'(20(2F10.6,2x))') (DIP_O(i,j),j=1,dim)
      end do
      write(u6,*)
      write(u6,'( 5x,a,i2,a,i3)') 'DIP_W, n = ',N,', m =',m
      write(u6,*)
      do i=1,dim
        write(u6,'(20(2F10.6,2x))') (DIP_W(i,j),j=1,dim)
      end do
    end if

    SP_DIPO(:) = czero
    SP_DIPW(:) = czero
    SP_MOW = trace(dim,DIP_O,DIP_W)
    do l=1,3
      SP_DIPO(l) = trace(dim,AMS(l,1:dim,1:dim),DIP_O(1:dim,1:dim))
      SP_DIPW(l) = trace(dim,AMS(l,1:dim,1:dim),DIP_W(1:dim,1:dim))

      B(l,n,-m) = SP_DIPO(l)/SP_MOW
      B(l,n,m) = SP_DIPW(l)/SP_MOW
    end do ! l
  end do ! m
end do ! n

BNMC(:,:,:) = cZero
BNMS(:,:,:) = cZero
do n=1,dim-1
  do m=0,N
    do l=1,3
      if (M == 0) then
        BNMC(l,n,m) = Half*(B(l,n,m)+B(l,n,-m))
      else
        m_fact = (-One)**M
        BNMC(l,n,m) = B(l,n,m)+m_fact*B(l,n,-m)
        BNMS(l,n,m) = -Onei*(B(l,n,m)-m_fact*B(l,n,-m))
      end if
    end do
  end do
end do !n

if (iprint > 2) then
  write(u6,'(100A)') ('-',i=1,47),'|'
  write(u6,'(A)') '  n  |  m  |   |       B       |       C       |'
  do N=1,dim-1
    write(u6,'(A)') '-----|-----|---|---------------|---------------|'
    do M=0,N
      if (M /= 0) write(u6,'(A)') '     |-----|---|---------------|---------------|'
      write(u6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))') N,'|',M,'|','X','|',real(BNMC(1,N,M)),'|',real(BNMS(1,N,M)),'|'
      write(u6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))') N,'|',M,'|','Y','|',real(BNMC(2,N,M)),'|',real(BNMS(2,N,M)),'|'
      write(u6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))') N,'|',M,'|','Z','|',real(BNMC(3,N,M)),'|',real(BNMS(3,N,M)),'|'
    end do
  end do
  write(u6,'(100A)') ('-',i=1,47),'|'
end if

HCF2(:,:,:,:) = cZero
do N=1,dim-1
  do M=0,N
    DIP_O(:,:) = cZero
    DIP_W(:,:) = cZero

    call Stewens_matrixel(N,M,dim,DIP_O,DIP_W,iprint)

    do l=1,3
      do i=1,dim
        do j=1,dim
          if (M == 0) then
            HCF2(N,l,i,j) = HCF2(N,l,i,j)+BNMC(l,N,M)*DIP_O(i,j)
          else
            m_fact = (-One)**M
            O1 = Half*(m_fact*DIP_W(i,j)+DIP_O(i,j))
            O2 = -Half*Onei*(m_fact*DIP_W(i,j)-DIP_O(i,j))

            HCF2(N,l,i,j) = HCF2(N,l,i,j)+BNMC(l,N,M)*O1+BNMS(l,N,M)*O2
          end if
        end do
      end do
    end do
  end do
end do ! n

if (iprint > 2) then
  do N=1,dim-1,2
    write(u6,*)
    write(u6,'( 5x,a,I2,a)') 'HCF2(',N,',l,i,j)'
    do l=1,3
      write(u6,*)
      write(u6,'(a,i3)') 'PROJECTION =',l
      do i=1,dim
        write(u6,'(20(2F12.8,2x))') (HCF2(N,l,i,j),j=1,dim)
      end do
    end do
    write(u6,*)
  end do
end if

gtens(:) = Zero
maxes(:,:) = Zero
call ATENS(HCF2(order,:,:,:),dim,gtens,maxes,1)

return

end subroutine mu_order
