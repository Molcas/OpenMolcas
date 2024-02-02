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

implicit none
integer, parameter :: wp = kind(0.d0)
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
real(kind=8) :: maxes(3,3)
complex(kind=8) :: DIP_O(dim,dim), DIP_W(dim,dim), B(3,dim,-dim:dim), BNMC(3,dim,0:dim), BNMS(3,dim,0:dim), SP_MOW, SP_DIPO(3), &
                   SP_DIPW(3), O1, O2, m_fact, trace
external :: trace
!------------------------------------------------------------

Z = (0.0_wp,0.0_wp)
! get the local pseudospin:
call pseudospin(MM,dim,Z,3,1,iprint)
! re-write MM and MS to the new pseudospin basis:
AMS = (0.0_wp,0.0_wp)
AMM = (0.0_wp,0.0_wp)
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
B(1:3,1:dim,-dim:dim) = (0.0_wp,0.0_wp)
do N=1,dim-1
  do M=0,N
    DIP_O = (0.0_wp,0.0_wp)
    DIP_W = (0.0_wp,0.0_wp)
    call Stewens_matrixel(N,M,dim,DIP_O,DIP_W,iprint)
    if (iprint > 5) then
      write(6,*)
      write(6,'( 5x,a,i2,a,i3)') 'DIP_O, n = ',N,', m =',m
      write(6,*)
      do i=1,dim
        write(6,'(20(2F10.6,2x))') (DIP_O(i,j),j=1,dim)
      end do
      write(6,*)
      write(6,'( 5x,a,i2,a,i3)') 'DIP_W, n = ',N,', m =',m
      write(6,*)
      do i=1,dim
        write(6,'(20(2F10.6,2x))') (DIP_W(i,j),j=1,dim)
      end do
    end if

    SP_DIPO = (0.0_wp,0.0_wp)
    SP_DIPW = (0.0_wp,0.0_wp)
    SP_MOW = (0.0_wp,0.0_wp)
    SP_MOW = trace(dim,DIP_O,DIP_W)
    do l=1,3
      SP_DIPO(l) = trace(dim,AMS(l,1:dim,1:dim),DIP_O(1:dim,1:dim))
      SP_DIPW(l) = trace(dim,AMS(l,1:dim,1:dim),DIP_W(1:dim,1:dim))

      B(l,n,-m) = SP_DIPO(l)/SP_MOW
      B(l,n,m) = SP_DIPW(l)/SP_MOW
    end do ! l
  end do ! m
end do ! n

BNMC = (0.0_wp,0.0_wp)
BNMS = (0.0_wp,0.0_wp)
do n=1,dim-1
  do m=0,N
    do l=1,3
      if (M == 0) then
        BNMC(l,n,m) = (0.5_wp,0.0_wp)*(B(l,n,m)+B(l,n,-m))
      else
        m_fact = cmplx((-1)**M,0,wp)
        BNMC(l,n,m) = B(l,n,m)+m_fact*B(l,n,-m)
        BNMS(l,n,m) = (B(l,n,m)-m_fact*B(l,n,-m))*(0.0_wp,-1.0_wp)
      end if
    end do
  end do
end do !n

if (iprint > 2) then
  write(6,'(100A)') ('-',i=1,47),'|'
  write(6,'(A)') '  n  |  m  |   |       B       |       C       |'
  do N=1,dim-1
    write(6,'(A)') '-----|-----|---|---------------|---------------|'
    do M=0,N
      if (M /= 0) write(6,'(A)') '     |-----|---|---------------|---------------|'
      write(6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))') N,'|',M,'|','X','|',dble(BNMC(1,N,M)),'|',dble(BNMS(1,N,M)),'|'
      write(6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))') N,'|',M,'|','Y','|',dble(BNMC(2,N,M)),'|',dble(BNMS(2,N,M)),'|'
      write(6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))') N,'|',M,'|','Z','|',dble(BNMC(3,N,M)),'|',dble(BNMS(3,N,M)),'|'
    end do
  end do
  write(6,'(100A)') ('-',i=1,47),'|'
end if

HCF2 = (0.0_wp,0.0_wp)
do N=1,dim-1
  do M=0,N
    DIP_O = (0.0_wp,0.0_wp)
    DIP_W = (0.0_wp,0.0_wp)

    call Stewens_matrixel(N,M,dim,DIP_O,DIP_W,iprint)

    do l=1,3
      do i=1,dim
        do j=1,dim
          if (M == 0) then
            HCF2(N,l,i,j) = HCF2(N,l,i,j)+BNMC(l,N,M)*DIP_O(i,j)
          else
            m_fact = (0.0_wp,0.0_wp)
            O1 = (0.0_wp,0.0_wp)
            O2 = (0.0_wp,0.0_wp)

            m_fact = cmplx((-1)**M,0,wp)
            O1 = (0.5_wp,0.0_wp)*(m_fact*DIP_W(i,j)+DIP_O(i,j))
            O2 = (0.0_wp,-0.5_wp)*(m_fact*DIP_W(i,j)-DIP_O(i,j))

            HCF2(N,l,i,j) = HCF2(N,l,i,j)+BNMC(l,N,M)*O1+BNMS(l,N,M)*O2
          end if
        end do
      end do
    end do
  end do
end do ! n

if (iprint > 2) then
  do N=1,dim-1,2
    write(6,*)
    write(6,'( 5x,a,I2,a)') 'HCF2(',N,',l,i,j)'
    do l=1,3
      write(6,*)
      write(6,'(a,i3)') 'PROJECTION =',l
      do i=1,dim
        write(6,'(20(2F12.8,2x))') (HCF2(N,l,i,j),j=1,dim)
      end do
    end do
    write(6,*)
  end do
end if

gtens = 0.0_wp
maxes = 0.0_wp
call ATENS(HCF2(order,:,:,:),dim,gtens,maxes,1)

return

end subroutine mu_order
