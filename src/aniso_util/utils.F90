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

subroutine prMom(a,m,n)

use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: a
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: M(3,n,n)
integer(kind=iwp) :: i, j, l
character(len=50) :: fmtline
character, parameter :: proj(3) = ['X','Y','Z']

write(u6,*)
write(u6,'(2a)') 'print: ',a
write(fmtline,'(a,i2,a)') '(',n,'(2f9.4,1x))'
do l=1,3
  write(u6,'(2a)') 'projection: ',proj(l)
  do i=1,n
    write(u6,fmtline) (M(l,i,j),j=1,n)
  end do
  write(u6,*)
end do

return

end subroutine prMom

subroutine prMom_herm(a,m,n)

use Constants, only: Three
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: a
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: M(3,n,n)
integer(kind=iwp) :: i, j, l
real(kind=wp) :: R

write(u6,*)
write(u6,'(2a)') 'print: ',a
do i=1,n
  do j=1,i
    R = (abs(M(1,i,j))+abs(M(2,i,j))+abs(M(3,i,j)))/Three
    write(u6,'(A,2I3,A,3(2F16.7,2x), 2F20.7)') 'i j: ',i,j,' <i|O|j>=',(M(l,i,j),l=1,3),R
  end do
  write(u6,*)
end do

return

end subroutine prMom_herm

subroutine pa_prMat(a,m,n)

use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: a
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: M(n,n)
integer(kind=iwp) :: i, j
character(len=50) fmtline

write(u6,*)
write(u6,'(2a)') 'print: ',a
write(fmtline,'(a,i2,a)') '(',n,'(2f12.4,1x))'
do i=1,n
  write(u6,fmtline) (M(i,j),j=1,n)
end do

return

end subroutine pa_prMat

subroutine pa_prMatR(a,m,n)

use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: a
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: M(n,n)
integer(kind=iwp) :: i, j
character(len=50) fmtline

write(u6,*)
write(u6,'(2a)') 'print: ',a
write(fmtline,'(a,i2,a)') '(',n,'(f19.14,1x))'
do i=1,n
  write(u6,fmtline) (M(i,j),j=1,n)
end do

return

end subroutine pa_prMatR

subroutine print_ZFS(a,m,n)

use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: a
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: M(n,n)
integer(kind=iwp) :: i, j, jEnd, k

write(u6,'(/)')
write(u6,'(A)') repeat('-',87)
write(u6,'(A)') a
do j=1,n,4
  jEnd = min(n,J+3)
  if (mod(n,2) == 0) then
    ! code for odd N
    write(u6,'(52A)') repeat('-',10),(repeat('-',24),k=j,jEnd),'|'
    write(u6,'(10x,A,50(8x,A,I3,A,7x,A))') '|',('|',2*i-n-1,'/2 >','|',i=j,jEnd)
    write(u6,'(52A)') repeat('-',10),'|',('---- Real ----- Imag --|',k=j,jEnd)
    ! print the matrix
    do i=1,n
      write(u6,'(1x,A,I3,A,1x,A,50(2F11.5,1x,A))') '<',2*i-n-1,'/2','| |',(M(k,i),'|',k=j,jEnd)
    end do
    write(u6,'(52A)') repeat('-',10),(repeat('-',24),k=j,jEnd),'|'

  else
    ! code for odd N

    write(u6,'(52A)') repeat('-',8),(repeat('-',24),k=j,jEnd),'|'
    write(u6,'(8x,A,50(8x,A,I3,A,9x,A))') '|',('|',-(n-1)/2-1+i,' >','|',i=j,jEnd)
    write(u6,'(52A)') repeat('-',8),'|',('---- Real ----- Imag --|',k=j,jEnd)
    do i=1,n
      write(u6,'(1x,A,I3,1x,A,50(2F11.5,1x,A))') '<',-(n-1)/2-1+i,'| |',(M(k,i),'|',k=j,jEnd)
    end do
    write(u6,'(52A)') repeat('-',8),(repeat('-',24),k=j,jEnd),'|'
  end if
end do !j

return

end subroutine print_ZFS

subroutine print_ZFS_eigenvectors(a,m,n)

use Definitions, only: wp, iwp, u6

implicit none
character, intent(in) :: a
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: M(n,n) ! eigenvectors
integer(kind=iwp) :: i, j, jEnd, k

write(u6,'(/)')
do j=1,n,2
  jEnd = min(n,j+1) ! '|  > |'
  write(u6,'(A,6A)') '--------|',('-----------------------------|',i=j,jEnd)
  write(u6,'(3A,6(6x,a,i3,5x,a))') ' | ',a,'M > |',('ab initio state',i,'|',i=j,jEnd)
  write(u6,'(A,6A)') '--------|',('-- Real ---- Imag --|-Weight-|',i=j,jEnd)

  do i=1,n
    if (mod(n,2) == 1) then
      write(u6,'(1x,A,1x,i2,A,6(2(ES22.14,1x),a,F6.1,1x,a))') '|',-(n-1)/2-(1-i),' > |',(real(M(i,k)),aimag(M(i,k)),'|', &
                                                              100.0_wp*(real(M(i,k))**2+aimag(M(i,k))**2),'%|',k=j,jEnd)
    else
      write(u6,'(A,i3,a,a,6(2(ES22.14,1x),a,F6.1,1x,a))') '|',-(n-1)+2*(i-1),'/2> ','|',(real(M(i,k)),aimag(M(i,k)),'|', &
                                                          100.0_wp*(real(M(i,k))**2+aimag(M(i,k))**2),'%|',k=j,jEnd)
    end if
  end do  !i

  write(u6,'(A,6A)') '--------|',('-----------------------------|',i=j,jEnd)
end do ! j

return

end subroutine print_ZFS_eigenvectors

subroutine print_ZFS_naoya(A,M,N)

use Definitions, only: wp, iwp, u6

implicit none
character, intent(in) :: a
integer(kind=iwp), intent(in) :: N ! dimension of the pseudospin
complex(kind=wp), intent(in) :: M(n,n) ! complex parameters to print
integer(kind=iwp) :: i, j, jEnd, k
real(kind=wp) :: Mi(n), Mr(n), Weight(n)

write(u6,'(/)')
do j=1,n,2
  jEnd = min(n,j+1) ! '|  > |'
  if (j == 1) write(u6,'(13A)') '--------|',(repeat('-',58),'|',i=j,jEnd)
  write(u6,'(3A,6(16x,a,i3,24x,a))') ' | ',a,'M > |',('ab initio state',i,'|',i=j,jEnd)
  write(u6,'(A,6A)') '--------|',('-------  Real  -------|------  Imaginary  -------|-Weight-|',i=j,jEnd)

  do i=1,n
    ! fix the sign:
    do k=j,jEnd
      Mr(k) = real(M(i,k))
      Mi(k) = aimag(M(i,k))
      Weight(k) = 100.0_wp*(Mr(k)*Mr(k)+Mi(k)*Mi(k))
    end do

    ! print it
    if (mod(n,2) == 1) then
      write(u6,'(1x,A,1x,i2,A,2(SP,2(1x,ES21.14,1x),a,S,F6.1,1x,a))') '|',-(n-1)/2-(1-i),' > |', &
                                                                      (Mr(k),Mi(k),'*I |',Weight(k),'%|',k=j,jEnd)
    else
      write(u6,'(A,i3,a,a,2(SP,2(1x,ES21.14,1x),a,S,F6.1,1x,a))') '|',-(n-1)+2*(i-1),'/2> ','|', &
                                                                 (Mr(k),Mi(k),'*I |',Weight(k),'%|',k=j,jEnd)
    end if
  end do ! i

  write(u6,'(13A)') '--------|',(repeat('-',58),'|',i=j,jEnd)
end do ! j

return

end subroutine print_ZFS_naoya

subroutine print_CFP_alpha(nlanth,n,B,C)

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nlanth, n  ! number of the lanthanide, dimension of the pseudospin
real(kind=wp), intent(in) :: B(n,0:n), C(n,0:n) ! real and imaginary CF parameters
!logical(kind=iwp), intent(in), optional :: print_all
integer(kind=iwp) :: k, q
real(kind=wp) :: a(6)

call set_an(nlanth,a)

!O(m1,m2) = Half*(Om+mQ*Op)
!W(m1,m2) = Half*Onei*(Om-mQ*Op)
write(u6,'(/)')
write(u6,'(A)') repeat('*',80)
write(u6,'(A)') 'The Crystal-Field Hamiltonian:'
write(u6,'(A)') '   Hcf = SUM_{k,q} alpha(k) * [ B(k,q) * O(k,q) +  C(k,q) * W(k,q) ];'
write(u6,'(A)') 'where:'
write(u6,'(A)') '   O(k,q) =  0.5 * ( (-1)**q * Y(k,+q) + Y(k,-q) );'
write(u6,'(A)') '   W(k,q) = -0.5 * ( (-1)**q * Y(k,+q) - Y(k,-q) ) * I;   (I = imaginary unit)'
write(u6,'(A)') '   k - the rank of the ITO, = 2, 4, 6;'
write(u6,'(A)') '   q - the component of the ITO, = 0, 1, 2, ... k;'
write(u6,'(A)') '   alpha(k) - Stevens coefficients;'
write(u6,'(A)') 'These operators have been defined in: '
write(u6,'(A)') '  L. F. Chibotaru, L.Ungur, J. Chem. Phys., 137, 064112 (2012).'

write(u6,'(2A)') repeat('-',76),'|'
write(u6,'(A)') '  k  |  q  |    1/alpha(k)  |         B(k,q)        |         C(k,q)        |'
do k=2,6,2
  if (abs(A(k)) > tiny(A)) then
    write(u6,'(A)') '-----|-----|----------------|-----------------------|-----------------------|'
    do q=0,k
      if (q == k/2) then
        write(u6,'(2(2x,I1,2x,A),F14.5,2x,A,2(ES22.14,1x,A))') k,'|',q,'|',One/a(k),'|',B(k,q)/a(k),'|',C(k,q)/a(k),'|'
      else
        write(u6,'(2(2x,I1,2x,A),16x,A,2(ES22.14,1x,A))') k,'|',q,'|','|',B(k,q)/a(k),'|',C(k,q)/a(k),'|'
      end if
    end do
  end if
end do
write(u6,'(2A)') repeat('-',76),'|'

return

end subroutine print_CFP_alpha

subroutine print_CFP_LCLU(n,B,C,print_all)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n ! dimension of the pseudospin
real(kind=wp), intent(in) :: B(n,0:n), C(n,0:n) ! real and imaginary CF parameters
logical(kind=iwp), intent(in) :: print_all
integer(kind=iwp) :: k, q

!O(m1,m2) = Half*(Om+mQ*Op)
!W(m1,m2) = Half*Onei*(Om-mQ*Op)
write(u6,'(/)')
write(u6,'(A)') repeat('*',80)
write(u6,'(A)') 'The Crystal-Field Hamiltonian:'
write(u6,'(A)') '   Hcf = SUM_{k,q} * [ B(k,q) * O(k,q) +  C(k,q) * W(k,q) ];'
write(u6,'(A)') 'where:'
write(u6,'(A)') '   O(k,q) =  0.5 * ( (-1)**q * Y(k,+q) + Y(k,-q) );'
write(u6,'(A)') '   W(k,q) = -0.5 * ( (-1)**q * Y(k,+q) - Y(k,-q) ) * I;   (I = imaginary unit)'
write(u6,'(A)') '   k - the rank of the ITO, = 2, 4, 6;'
write(u6,'(A)') '   q - the component of the ITO, = 0, 1, 2, ... k;'
write(u6,'(A)') 'These operators have been defined in: '
write(u6,'(A)') '  L. F. Chibotaru, L.Ungur, J. Chem. Phys., 137, 064112 (2012).'

write(u6,'(2A)') repeat('-',59),'|'
write(u6,'(A)') '  k  |  q  |         B(k,q)        |         C(k,q)        |'
if (print_all) then
  do k=2,n-1
    write(u6,'(A)') '-----|-----|-----------------------|-----------------------|'
    do q=0,k
      write(u6,'(2(1x,I2,2x,A),2(ES22.14,1x,A))') k,'|',q,'|',B(k,q),'|',C(k,q),'|'
    end do
  end do
else
  do k=2,n-1,2
    write(u6,'(A)') '-----|-----|-----------------------|-----------------------|'
    do q=0,k
      write(u6,'(2(1x,I2,2x,A),2(ES22.14,1x,A))') k,'|',q,'|',B(k,q),'|',C(k,q),'|'
    end do
  end do
end if
write(u6,'(2A)') repeat('-',59),'|'

return

end subroutine print_CFP_LCLU

subroutine print_CFP_stev(n,B,print_all)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n ! dimension of the pseudospin
real(kind=wp), intent(in) :: B(n,-n:n) ! real and imaginary CF parameters
logical(kind=iwp), intent(in) :: print_all
integer(kind=iwp) :: iq, k, kmax, q
real(kind=wp) :: f, knm(12,0:12)

call set_knm(knm)
write(u6,'(/)')
write(u6,'(A)') repeat('*',80)
write(u6,'(A)') 'The Crystal-Field Hamiltonian:'
write(u6,'(A)') '   Hcf = SUM_{k,q} * [ B(k,q) * O(k,q) ];'
write(u6,'(A)') 'where:'
write(u6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO) as defined in:'
write(u6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State Phys.,18(1985) 1415-1430.'
write(u6,'(10x,A)') '2. Implemented in the "EasySpin" function in MATLAB, www.easyspin.org.'
write(u6,'(A)') '   k - the rank of the ITO, = 2, 4, 6, 8, 10, 12.'
write(u6,'(A)') '   q - the component of the ITO, = -k, -k+1, ... 0, 1, ... k;'
if ((n-1) > 12) then
  write(u6,'(A)') 'k = 12 may not be the highest rank of the ITO for this case, but it '
  write(u6,'(A)') 'is the maximal k implemented in the "EasySpin" function in MATLAB.'
end if
write(u6,'(A)') 'Knm are proportionality coefficients between the ESO and operators defined in '
write(u6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
write(u6,'(2A)') repeat('-',48),'|'
write(u6,'(A)') '  k |  q  |    (K)^2    |         B(k,q)        |'
if ((n-1) > 12) then
  kmax = 12
else
  kmax = n-1
end if

if (print_all) then
  do k=2,kmax
    write(u6,'(A)') '----|-----|-------------|-----------------------|'
    do q=-k,k
      iq = abs(q)
      f = knm(k,iq)*knm(k,iq)
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,2(ES22.14,1x,A))') k,'|',q,'|',f,'|',B(k,q),'|'
    end do
  end do
else
  do k=2,kmax,2
    write(u6,'(A)') '----|-----|-------------|-----------------------|'
    do q=-k,k
      iq = abs(q)
      f = knm(k,iq)*knm(k,iq)
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,2(ES22.14,1x,A))') k,'|',q,'|',f,'|',B(k,q),'|'
    end do
  end do
end if
write(u6,'(2A)') repeat('-',48),'|'

return

end subroutine print_CFP_stev

subroutine print_CFP_naoya(N,A,print_all)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N ! dimension of the pseudospin
complex(kind=wp), intent(in) :: A(n-1,-(n-1):n-1) ! complex parameters to print
logical(kind=iwp), intent(in) :: print_all
integer(kind=iwp) :: k, q
real(kind=wp) :: Ai, Ar

write(u6,'(/)')
write(u6,'(A)') repeat('*',80)
write(u6,'(A)') 'The Crystal-Field Hamiltonian:'
write(u6,'(A)')
write(u6,'(A)') '   Hcf = SUM_{k,q} * [ B(k,q) * O(k,q) ];'
write(u6,'(A)')
write(u6,'(A)') ' where:                                  '
write(u6,'(A)') '   O(k,q) =  Irreducible Tensor Operators'
write(u6,'(A)') '             defined as follows:         '
write(u6,'(A)')
write(u6,'(A)') '          Y(k,q)             CG(J,M2,k,q,J,M1)'
write(u6,'(A)') ' < J,M1 | ------ | J,M2 >  = -----------------'
write(u6,'(A)') '          Y(k,0)              CG(J,J,k,0,J,J) '
write(u6,'(A)')
write(u6,'(A)') '  CG - Clebsh-Gordan Coefficient:'
write(u6,'(A)') '                            c,gm     '
write(u6,'(A)') '      CG(a,al,b,bt,c,gm) = C         '
write(u6,'(A)') '                            a,al,b,bt'
write(u6,'(A)')
write(u6,'(A)') '   k - the rank of the ITO, = 2, 4, 6, 8, 10, 12.'
write(u6,'(A)') '   q - the component of the ITO, = -k, -k+1, ... 0, 1, ... k;'
! thee table:
write(u6,'(2A)') repeat('-',59),'|'
write(u6,'(A,11x,A,12x,A)') '  k |  q  |','Complex parameter  A(k,q)','|'
write(u6,'(A)') '----|-----|--------  Real  ------|-----  Imaginary  -------|'

if (print_all) then
  do k=2,N-1
    do q=-k,k
      Ar = real(A(k,q))
      Ai = aimag(A(k,q))

      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),SP,ES21.14,1x,ES21.14,A)') k,'|',q,'| ',Ar,Ai,' *I |'

    end do
    if (k /= (N-1-mod(N-1,2))) write(u6,'(A)') '----|-----|----------------------|-------------------------|'
  end do
else
  do k=2,N-1,2
    do q=-k,k
      Ar = real(A(k,q))
      Ai = aimag(A(k,q))

      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),SP,ES21.14,1x,ES21.14,A)') k,'|',q,'| ', Ar,Ai,' *I |'

    end do
    if (k /= (N-1-mod(N-1,2))) write(u6,'(A)') '----|-----|----------------------|-------------------------|'
  end do
end if

write(u6,'(2A)') repeat('-',59),'|'

return

end subroutine print_cfp_naoya

subroutine print_MOM_ITO_stev(n,B,print_all)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n ! dimension of the pseudospin
real(kind=wp), intent(in) :: B(3,n,-n:n) ! real and imaginary CF parameters, for x,y,z
logical(kind=iwp), intent(in) :: print_all
integer(kind=iwp) :: iq, k, kmax, q
real(kind=wp) :: f, knm(12,0:12)

call set_knm(knm)
write(u6,'(/)')
write(u6,'(A)') repeat('*',80)
write(u6,'(A)') 'The magnetic moment is decomposed in Stev ITO:'
write(u6,'(A)') '   Hcf = SUM_{k,q} * [ B(k,q) * O(k,q) ];'
write(u6,'(A)') 'where:'
write(u6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO) as defined in:'
write(u6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State Phys.,18(1985) 1415-1430.'
write(u6,'(10x,A)') '2. Implemented in the "EasySpin" function in MATLAB, www.easyspin.org.'
write(u6,'(A)') '   k - the rank of the ITO, = 1, 3, 5, 7, 9, 11.'
write(u6,'(A)') '   q - the component of the ITO, = -k, -k+1, ... 0, 1, ... k;'
if ((n-1) > 12) then
  write(u6,'(A)') 'k = 12 may not be the highest rank of the ITO for this case, but it '
  write(u6,'(A)') 'is the maximal k implemented in the "EasySpin" function in MATLAB.'
end if
write(u6,'(A)') 'Knm are proportionality coefficients between the ESO and operators defined in '
write(u6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
write(u6,'(2A)') repeat('-',96),'|'
write(u6,'(A)') '  k |  q  |    (K)^2    |        B(k,q) - X     |        B(k,q) - Y     |        B(k,q) - Z     |'
if ((n-1) > 12) then
  kmax = 12
else
  kmax = n-1
end if

if (print_all) then

  do k=1,kmax
    write(u6,'(A)') '----|-----|-------------|-----------------------|-----------------------|-----------------------|'
    do q=-k,k
      iq = abs(q)
      f = knm(k,iq)*knm(k,iq)
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,3(ES22.14,1x,A))') k,'|',q,'|',f,'|',B(1,k,q),'|',B(2,k,q),'|',B(3,k,q),'|'
    end do
  end do

else

  do k=1,kmax,2
    write(u6,'(A)') '----|-----|-------------|-----------------------|-----------------------|-----------------------|'
    do q=-k,k
      iq = abs(q)
      f = knm(k,iq)*knm(k,iq)
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,3(ES22.14,1x,A))') k,'|',q,'|',f,'|',B(1,k,q),'|',B(2,k,q),'|',B(3,k,q),'|'
    end do
  end do

end if
write(u6,'(2A)') repeat('-',96),'|'

return

end subroutine print_MOM_ITO_stev
