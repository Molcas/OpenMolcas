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

subroutine ITO(N,k,q,C0,Cp,Cm)
! This soubroutine implements the core ITO generation
! using the following ITO convention:
!                               S,M2
!                             CG
!                               S,M1,k,q
!  < S,M1 | O_k_q | S,M2 > = -------------
!                               S,S
!                             CG
!                               S,S,k,0
!
! The denominator is computed using a simpler formula:
! Varshalovich, p252, formula (42):
!
! N = 2S+1, dimension   (input)
! k = rank of ITO       (input)
! q = projection of ITO (input)
! Cp = O+ operator (output), complex
! Cm = O- operator (output), complex
! C0 = CG0  (output), real number, positive

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: n, k, q
real(kind=8), intent(out) :: C0
complex(kind=8), intent(out) :: Cp(n,n), Cm(n,n)
! local
integer :: m1, m2
real(kind=8) :: rm1, rm2, rS, rK, rQ, CGp, CGm, fct
external :: fct

call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cp,1)
call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cm,1)

rS = dble(n-1)/2.0_wp
rK = dble(k)
rQ = dble(q)
C0 = fct(n-1)*sqrt(dble(n)/(fct(n-k-1)*fct(n+k)))
! M1 and M2 go from max value to min value
do m1=1,n
  do m2=1,n
    rm1 = rS-dble(m1-1)
    rm2 = rS-dble(m2-1)
    call Clebsh_Gordan(rS,rm2,rK,rQ,rS,rm1,CGp)
    call Clebsh_Gordan(rS,rm2,rK,-rQ,rS,rm1,CGm)
    Cp(m1,m2) = cmplx(CGp/C0,0.0_wp,wp)
    Cm(m1,m2) = cmplx(CGm/C0,0.0_wp,wp)
  end do
end do

return

end subroutine ITO
!=!=
subroutine Liviu_ITO(n,k,q,O,W,redME)

! generate O, W as in previous Stewens_matrixel function:
! redME =  ratio between Naoya's and Liviu operators
implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: n, k, q
complex(kind=8), intent(out) :: O(n,n), W(n,n), redME
! local
real(kind=8) :: CR, C0

call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,O,1)
call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,W,1)
CR = 0.0_wp
C0 = 0.0_wp
redME = (0.0_wp,0.0_wp)
call coeff_redus_sub(n,k,CR)
call ITO(n,k,q,C0,O,W)
redME = cmplx(C0*CR,0.0_wp,wp)
call zscal_(n*n,redME,W,1)
call zscal_(n*n,redME,O,1)

return

end subroutine Liviu_ITO
!=!=
subroutine Stev_ITO(n,k,q,O,W,redME)
! generate O, W as in previous Stewens_matrixel function:
! redME =  ratio between Naoya's and Liviu operators
! scaled as for Stevens Operators by knm

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: n, k, q
complex(kind=8), intent(out) :: O(n,n), W(n,n), redME
! local
real(kind=8) :: F, knm(12,0:12)

call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,O,1)
call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,W,1)
redME = (0.0_wp,0.0_wp)
if ((k > 12) .or. (q > 12)) return

call set_knm(knm)
call Liviu_ITO(n,k,q,O,W,redME)
F = 1.0_wp/knm(k,q)
redME = cmplx(F,0.0_wp,wp)
call zscal_(n*n,redME,W,1)
call zscal_(n*n,redME,O,1)

return

end subroutine Stev_ITO
!=!=
subroutine ESO(n,k,q,O,W,redME)
! generate Hermitian ESO operators, as in MATLAB EasySpin(stev) function
! only for q >= 0

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: n, k, q
complex(kind=8), intent(out) :: O(n,n), W(n,n), redME
! local
integer :: m1, m2
real(kind=8) :: CR, C0, knm(12,0:12), F
complex(kind=8) :: mQ, HALF_R, FALF_I
complex(kind=8), allocatable :: Cp(:,:), Cm(:,:)

call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,O,1)
call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,W,1)
redME = (0.0_wp,0.0_wp)
if ((k > 12) .or. (q > 12)) return ! not available in MATLAB

call mma_allocate(Cp,n,n,'Cp')
call mma_allocate(Cm,n,n,'Cm')
call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cp,1)
call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cm,1)

call set_knm(knm)
call coeff_redus_sub(n,k,CR)
call ITO(n,k,q,C0,Cp,Cm)

F = C0*CR/knm(k,q)
redME = cmplx(F,0.0_wp,wp)
mQ = cmplx((-1)**q,0.0_wp,wp)
HALF_R = (0.5_wp,0.0_wp)
FALF_I = (0.0_wp,0.5_wp)
do m1=1,n
  do m2=1,n
    O(m1,m2) = HALF_R*redME*(Cm(m1,m2)+mQ*Cp(m1,m2))
    W(m1,m2) = FALF_I*redME*(Cm(m1,m2)-mQ*Cp(m1,m2))
  end do
end do
call mma_deallocate(Cp)
call mma_deallocate(Cm)

return

end subroutine ESO
!=!=
subroutine Liviu_ESO(n,k,q,O,W,redME)
! generate Hermitian ESO operators, as in MATLAB EasySpin(stev) function
! only for q >= 0
! Om = Ok-q
! Op = Ok+q
! the only difference with MATLAB's ESO are scaling factor Knm

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: n, k, q
complex(kind=8), intent(out) :: O(n,n), W(n,n), redME
! local
integer :: m1, m2
real(kind=8) :: CR, C0
complex(kind=8) :: Om, Op, mQ
complex(kind=8), allocatable :: Cp(:,:), Cm(:,:)

call mma_allocate(Cp,n,n,'Cp')
call mma_allocate(Cm,n,n,'Cm')
call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cp,1)
call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cm,1)

call coeff_redus_sub(n,k,CR)
call ITO(n,k,q,C0,Cp,Cm)
redME = cmplx(C0*CR,0.0_wp,wp)
mQ = cmplx(dble((-1)**q),0.0_wp,wp)
do m1=1,n
  do m2=1,n
    Om = Cm(m1,m2)*redME
    Op = Cp(m1,m2)*redME
    O(m1,m2) = (0.5_wp,0.0_wp)*(Om+mQ*Op)
    W(m1,m2) = (0.0_wp,0.5_wp)*(Om-mQ*Op)
  end do
end do
call mma_deallocate(Cp)
call mma_deallocate(Cm)

return

end subroutine Liviu_ESO
!=!=
subroutine Stewens_matrixel(N,M,dim,ITO_O,ITO_W,IPRINT)
! This Subroutine calculates the matrix elements of the ITO  (On)
! on the basis of the eigenfunctions of the effective spin.
!
! N -- the rank of the ITO (On)
! dim -- the multiplicity of the effective spin
! Dip_Stewens(N, L, dim, dim) -- the matrix elements of the ITO tensor
!               operators in the basis of effective spin eigenfunctions

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: N, M, dim, IPRINT
complex(kind=8), intent(out) :: ITO_O(dim,dim), ITO_W(dim,dim)
integer :: npar, i, j, ms1, ms2
real(kind=8) :: a, al, b, bt, c, gm, COEFF_REDUS, coeffCG
complex(kind=8) :: ITO_PLUS(-dim:dim,-dim:dim), ITO_MINUS(-dim:dim,-dim:dim)
!***********************************************************************

NPAR = mod(dim,2)
COEFF_REDUS = 0.0_wp
ITO_PLUS(-dim:dim,-dim:dim) = (0.0_wp,0.0_wp)
ITO_MINUS(-dim:dim,-dim:dim) = (0.0_wp,0.0_wp)
ITO_O(1:dim,1:dim) = (0.0_wp,0.0_wp)
ITO_W(1:dim,1:dim) = (0.0_wp,0.0_wp)

call COEFF_REDUS_SUB(dim,N,COEFF_REDUS)
!COEFF_REDUS = 1.0_wp
a = dble(N)
al = dble(M)
c = (dble(dim)-1.0_wp)/2.0_wp
b = (dble(dim)-1.0_wp)/2.0_wp
do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
  if ((ms1 == 0) .and. (NPAR == 0)) go to 160
  if (NPAR == 0) then
    if (ms1 < 0) then
      gm = dble(ms1)+0.5_wp
    else
      gm = dble(ms1)-0.5_wp
    end if
  else
    gm = dble(ms1)
  end if

  do ms2=-(dim-NPAR)/2,(dim-NPAR)/2
    if ((ms2 == 0) .and. (NPAR == 0)) go to 150
    if (NPAR == 0) then
      if (ms2 < 0) then
        bt = dble(ms2)+0.5_wp
      else
        bt = dble(ms2)-0.5_wp
      end if
    else
      bt = dble(ms2)
    end if
    coeffCG = 0.0_wp

    call Clebsh_Gordan(a,al,b,bt,c,gm,coeffCG)

    ITO_PLUS(ms1,ms2) = cmplx(coeffCG*COEFF_REDUS,0.0_wp,wp)

    if (iprint > 5) &
      write(6,'(5x,2(a,i3,2x),6(a,f4.1,2x),2(a,f14.10,2x))') 'ms1=',ms1,'ms2=',ms2,'a=',a,'al=',al,'b=',b,'bt=',bt,'c=',c, &
                                                             'gm=',gm,'coeffCG=',coeffCG,'coeffCG^2=',coeffCG**2
150 continue
  end do ! ms2
160 continue
end do ! ms1

al = -dble(M)

do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
  if ((ms1 == 0) .and. (NPAR == 0)) go to 161
  if (NPAR == 0) then
    if (ms1 < 0) then
      gm = dble(ms1)+0.5_wp
    else
      gm = dble(ms1)-0.5_wp
    end if
  else
    gm = dble(ms1)
  end if

  do ms2=-(dim-NPAR)/2,(dim-NPAR)/2
    if ((ms2 == 0) .and. (NPAR == 0)) go to 151

    if (NPAR == 0) then
      if (ms2 < 0) then
        bt = dble(ms2)+0.5_wp
      else
        bt = dble(ms2)-0.5_wp
      end if
    else
      bt = dble(ms2)
    end if
    coeffCG = 0.0_wp

    call Clebsh_Gordan(a,al,b,bt,c,gm,coeffCG)

    ITO_MINUS(ms1,ms2) = cmplx(coeffCG*COEFF_REDUS,0.0_wp,wp)

    if (iprint > 5) &
      write(6,'(5x,2(a,i3,2x),6(a,f4.1,2x),2(a,f14.10,2x))') 'ms1=',ms1,'ms2=',ms2,'a=',a,'al=',al,'b=',b,'bt=',bt,'c=',c, &
                                                             'gm=',gm,'coeffCG=',coeffCG,'coeffCG^2=',coeffCG**2
151 continue
  end do
161 continue
end do

i = 0
do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
  if ((ms1 == 0) .and. (NPAR == 0)) go to 180
  i = i+1
  j = 0

  do ms2=-(dim-NPAR)/2,(dim-NPAR)/2
    if ((ms2 == 0) .and. (NPAR == 0)) go to 170
    j = j+1
    ITO_O(i,j) = ITO_PLUS(ms1,ms2)
    ITO_W(i,j) = ITO_MINUS(ms1,ms2)
170 continue
  end do
180 continue
end do

if (IPRINT > 3) then
  write(6,'(/)')
  write(6,'(5X,A,2I3)') 'Operator ITO_PLUS',N,M
  write(6,*)
  do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
    if ((ms1 == 0) .and. (NPAR == 0)) go to 158
    if (NPAR == 1) write(6,'(16(2X,2ES12.3))') (ITO_PLUS(ms1,ms2),ms2=-(dim-NPAR)/2,(dim-NPAR)/2)
    if (NPAR == 0) write(6,'(16(2X,2ES12.3))') (ITO_PLUS(ms1,ms2),ms2=-dim/2,-1),(ITO_PLUS(ms1,ms2),ms2=1,dim/2)
158 continue
  end do ! ms1

  write(6,'(/)')
  write(6,'(5X,A,2I3)') 'Operator ITO_MINUS',N,M
  write(6,*)
  do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
    if ((ms1 == 0) .and. (NPAR == 0)) go to 159
    if (NPAR == 1) write(6,'(16(2X,2ES12.3))') (ITO_MINUS(ms1,ms2),ms2=-(dim-NPAR)/2,(dim-NPAR)/2)
    if (NPAR == 0) write(6,'(16(2X,2ES12.3))') (ITO_MINUS(ms1,ms2),ms2=-dim/2,-1),(ITO_MINUS(ms1,ms2),ms2=1,dim/2)
159 continue
  end do ! ms1

  write(6,'(/////)')
  write(6,'(5X,A,2I3)') 'ITO_O',N,M
  write(6,*)
  do i=1,dim
    write(6,'(16(2X,2ES12.3))') (ITO_O(i,j),j=1,dim)
  end do   ! i
  write(6,'(/)')
  write(6,'(5X,A,2I3)') 'ITO_W',N,M
  write(6,*)
  do i=1,dim
    write(6,'(16(2X,2ES12.3))') (ITO_W(i,j),j=1,dim)
  end do ! i
end if ! iPrint

return

end subroutine Stewens_matrixel
!=!=
real*8 function fct(n)
! this function provides correct answer till n=169 only

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: n
integer :: i
real(kind=8) :: xct

xct = 1.0_wp
fct = 1.0_wp
if (n < 0) then
  write(6,'(A,i0)') 'FCT:  N<0 !'
  write(6,'(A,i0)') 'N = ',N
  write(6,'(A   )') 'It is an impossible case.'
  fct = -9.d99
  return

else if (n == 0) then
  return

else if ((n > 0) .and. (n <= 169)) then
  do i=1,n
    xct = xct*dble(i)
  end do

else
  write(6,'(A,i0)') 'FCT:   N = ',N
  write(6,'(A)') 'Factorial of N>169 overflows on x86_64'
  write(6,'(A)') 'Use higher numerical precision, or rethink your algorithm.'
end if

fct = xct

return

end function fct
!=!=
subroutine COEFF_REDUS_SUB(dim,N,COEFF_REDUS)
! THIS Subroutine ReturnS THE VALUE OF THE REDUCED MATRIX ELEMENT  <S|| On ||S>
! WHERE S-EFFECTIVE SPIN, On -- THE HIGHER ORDER SPIN OPERATORS
!
! The convention is to follow the paper:
! C. Rudowicz and C.Y. Chung
! J. Phys.: Condens. Matter, 2004, 16, pp. 5825
! doi:10.1088/0953-8984/16/32/018
!
! with a minor addition of the Norm(N) factor which depends on the operator rank
! This norm makes the matrix elements of ESO identical to those of the
! operators from EasySpin (stev).

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: N, dim
real(kind=8), intent(out) :: COEFF_REDUS
integer :: i
real(kind=8) :: FCT, Norm(100), s1, s2
external :: fct

COEFF_REDUS = 0.0_wp
do i=1,100
  Norm(i) = 0.0_wp
end do
Norm(1) = 1.0_wp
Norm(2) = 2.0_wp
Norm(3) = 2.0_wp
Norm(4) = 8.0_wp
Norm(5) = 8.0_wp
Norm(6) = 16.0_wp
Norm(7) = 16.0_wp
Norm(8) = 128.0_wp
Norm(9) = 128.0_wp
Norm(10) = 256.0_wp
Norm(11) = 256.0_wp
Norm(12) = 1024.0_wp
Norm(13) = 1024.0_wp
Norm(14) = 2048.0_wp
Norm(15) = 2048.0_wp
Norm(16) = 32768.0_wp
Norm(17) = 32768.0_wp
Norm(18) = 65536.0_wp
Norm(19) = 65536.0_wp
Norm(20) = 262144.0_wp
Norm(21) = 262144.0_wp
Norm(22) = 524288.0_wp
Norm(23) = 524288.0_wp
Norm(24) = 4194304.0_wp
Norm(25) = 4194304.0_wp
Norm(26) = 8388608.0_wp
Norm(27) = 8388608.0_wp
Norm(28) = 33554432.0_wp
Norm(29) = 33554432.0_wp
Norm(30) = 67108867.0_wp
Norm(31) = 67108867.0_wp
Norm(32) = 2147483648.0_wp
Norm(33) = 2147483648.0_wp
Norm(34) = 4294967296.0_wp
Norm(35) = 4294967296.0_wp
Norm(36) = 17179869184.0_wp
Norm(37) = 17179869184.0_wp
Norm(38) = 34359738368.0_wp
Norm(39) = 34359738368.0_wp
Norm(40) = 274877906944.0_wp
Norm(41) = 274877906944.0_wp
Norm(42) = 549755813888.0_wp
Norm(43) = 549755813888.0_wp
Norm(44) = 2199023255552.0_wp
Norm(45) = 2199023255552.0_wp
Norm(46) = 4398046511104.0_wp
Norm(47) = 4398046511104.0_wp
Norm(48) = 70368744177664.0_wp
Norm(49) = 70368744177664.0_wp
Norm(50) = 140737488355328.0_wp
Norm(51) = 140737488355328.0_wp
Norm(52) = 562949953421312.0_wp
Norm(53) = 562949953421312.0_wp
Norm(54) = 1125899906842624.0_wp
Norm(55) = 1125899906842624.0_wp
Norm(56) = 9007199254740992.0_wp
Norm(57) = 9007199254740992.0_wp
Norm(58) = 18014398509481984.0_wp
Norm(59) = 18014398509481984.0_wp
Norm(60) = 72057594037927936.0_wp
Norm(61) = 72057594037927936.0_wp
Norm(62) = 144115188075855872.0_wp
Norm(63) = 144115188075855872.0_wp
Norm(64) = 9223372036854775808.0_wp
Norm(65) = 9223372036854775808.0_wp
Norm(66) = 18446744073709551616.0_wp
Norm(67) = 18446744073709551616.0_wp
Norm(68) = 73786976294838206464.0_wp
Norm(69) = 73786976294838206464.0_wp
Norm(70) = 147573952589676412928.0_wp
Norm(71) = 147573952589676412928.0_wp
Norm(72) = 1180591620717411303424.0_wp
Norm(73) = 1180591620717411303424.0_wp
Norm(74) = 2361183241434822606848.0_wp
Norm(75) = 2361183241434822606848.0_wp
Norm(76) = 9444732965739290427392.0_wp
Norm(77) = 9444732965739290427392.0_wp
Norm(78) = 18889465931478580854784.0_wp
Norm(79) = 18889465931478580854784.0_wp
Norm(80) = 302231454903657293676544.0_wp
Norm(81) = 302231454903657293676544.0_wp
!  new method
s1 = 0.0_wp
s2 = 0.0_wp
s1 = sqrt(fct(dim+N)/fct(dim-N-1))
s2 = dble(2**N)*sqrt(dble(dim))

COEFF_REDUS = Norm(N)*s1/s2

return

end subroutine COEFF_REDUS_SUB
!=!=
subroutine Clebsh_Gordan(a,al,b,bt,c,gm,coeffCG)

implicit none
integer, parameter :: wp = kind(0.d0)
real(kind=8), intent(in) :: a, al, b, bt, c, gm
real(kind=8), intent(out) :: coeffCG
real(kind=8) :: u, fct, s1, s2
integer :: lb1, lb2, i
external :: fct

! exclude the cases for which CG coefficients are exactly zero
coeffCG = 0.0_wp
if ((al+bt) /= gm) return
if (a < 0.0_wp) return
if (b < 0.0_wp) return
if (c < 0.0_wp) return
if (abs(al) > a) return
if (abs(bt) > b) return
if (abs(gm) > c) return
if ((abs(a-b) > c) .or. ((a+b) < c)) return
if ((abs(b-c) > a) .or. ((b+c) < a)) return
if ((abs(c-a) > b) .or. ((c+a) < b)) return
if (mod(nint(2.0_wp*a),2) /= mod(nint(2.0_wp*abs(al)),2)) return
if (mod(nint(2.0_wp*b),2) /= mod(nint(2.0_wp*abs(bt)),2)) return
if (mod(nint(2.0_wp*c),2) /= mod(nint(2.0_wp*abs(gm)),2)) return
u = 0.0_wp
lb1 = int(min(c-b+al,c-a-bt))
lb2 = int(min(a+b-c,a-al,b+bt))
if (lb1 < 0) then
  if (-lb1 > lb2) then
    coeffCG = 0.0_wp
    return
  else

    do i=-lb1,lb2
      u = u+dble((-1)**i)/ &
          dble(fct(i)*fct(nint(a-al-i))*fct(nint(b+bt-i))*fct(nint(a+b-c-i))*fct(nint(c-b+al+i))*fct(nint(c-a-bt+i)))
    end do
  end if
else
  do i=0,lb2
    u = u+dble((-1)**i)/dble(fct(i)*fct(nint(a-al-i))*fct(nint(b+bt-i))*fct(nint(a+b-c-i))*fct(nint(c-b+al+i))*fct(nint(c-a-bt+i)))
  end do
end if

s1 = 0.0_wp
s2 = 0.0_wp
s1 = sqrt(dble(fct(nint(a+b-c))*fct(nint(a-b+c))*fct(nint(-a+b+c)))/dble(fct(nint(a+b+c+1))))

s2 = sqrt(dble(fct(nint(a+al))*fct(nint(a-al))*fct(nint(b+bt))*fct(nint(b-bt))*fct(nint(c+gm))*fct(nint(c-gm))*(2*c+1)))

coeffCG = u*s1*s2

return

end subroutine Clebsh_Gordan
!=!=
real*8 function W9J(a,b,c,d,e,f,g,h,j)

! Calculates a Wigner 9-j symbol. Argument a-j are Integer and are
! twice the true value of the 9-j's arguments, in the form
! { a b c }
! { d e f }
! { g h j }
! this is the implementation the formula 10.2.4. (20) from:
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: a, b, c, d, e, f, g, h, j
integer :: n, nlow, nhig
real(kind=8) :: W6J
logical :: check_triangle
external :: W6J, check_triangle

W9j = 0.0_wp
if (mod(a+b,2) /= mod(c,2)) return
if (mod(d+e,2) /= mod(f,2)) return
if (mod(g+h,2) /= mod(j,2)) return
if (mod(a+d,2) /= mod(g,2)) return
if (mod(b+e,2) /= mod(h,2)) return
if (mod(c+f,2) /= mod(j,2)) return
if ((abs(a-b) > c) .or. (a+b < c)) return
if ((abs(d-e) > f) .or. (d+e < f)) return
if ((abs(g-h) > j) .or. (g+h < j)) return
if ((abs(a-d) > g) .or. (a+d < g)) return
if ((abs(b-e) > h) .or. (b+e < h)) return
if ((abs(c-f) > j) .or. (c+f < j)) return
if (.not. check_triangle(a,b,c)) return
if (.not. check_triangle(d,e,f)) return
if (.not. check_triangle(g,h,j)) return
if (.not. check_triangle(a,d,g)) return
if (.not. check_triangle(b,e,h)) return
if (.not. check_triangle(c,f,j)) return

nlow = max(abs(a-j)/2,abs(d-h)/2,abs(b-f)/2)
nhig = min((a+j)/2,(d+h)/2,(b+f)/2)

do n=nlow,nhig
  W9j = W9j+dble(2*n+1)*dble((-1)**(2*n))*W6J(a,b,c,f,j,2*n)*W6J(d,e,f,b,2*n,h)*W6J(g,h,j,2*n,a,d)
end do

return

end function W9J
!=!=
real*8 function W6J(a,b,c,d,e,f)
! Calculates a Wigner 6-j symbol. Argument a-f are positive Integer
! and are twice the true value of the 6-j's arguments, in the form
! { a b c }
! { d e f }
!
! this is the implementation the formula 9.2.1. (1) from:
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: a, b, c, d, e, f
integer :: n, nlow, nhig
real(kind=8) :: dlt, sum, fct, isum
logical :: check_triangle
external :: fct, dlt, check_triangle

W6J = 0.0_wp
if (mod(a+b,2) /= mod(c,2)) return
if (mod(c+d,2) /= mod(e,2)) return
if (mod(a+e,2) /= mod(f,2)) return
if (mod(b+d,2) /= mod(f,2)) return
if ((abs(a-b) > c) .or. (a+b < c)) return
if ((abs(c-d) > e) .or. (c+d < e)) return
if ((abs(a-e) > f) .or. (a+e < f)) return
if ((abs(b-d) > f) .or. (b+d < f)) return

if (.not. check_triangle(a,b,c)) return
if (.not. check_triangle(c,d,e)) return
if (.not. check_triangle(a,e,f)) return
if (.not. check_triangle(b,d,f)) return

nlow = 0
nhig = 0
nlow = max((a+b+c)/2,(c+d+e)/2,(b+d+f)/2,(a+e+f)/2)
nhig = min((a+b+d+e)/2,(b+c+e+f)/2,(a+c+d+f)/2)

sum = 0.0_wp
do n=nlow,nhig
  isum = dble((-1)**n)*fct(n+1)/fct((a+c+d+f)/2-n)/fct((b+c+e+f)/2-n)/fct(n-(a+b+c)/2)/fct(n-(c+d+e)/2)/fct(n-(a+e+f)/2)/ &
         fct(n-(b+d+f)/2)/fct((a+b+d+e)/2-n)
  sum = sum+isum
end do
W6J = dlt(a,b,c)*dlt(c,d,e)*dlt(a,e,f)*dlt(b,d,f)*sum

return

end function W6J
!=!=
real*8 function W3J(j1,j2,j3,m1,m2,m3)
! Calculates a Wigner 3-j symbol. Argument j1,j2,j3 are positive Integer
! and are twice the true value of the 3-j's arguments, in the form
! { j1 j2 j3 }
! { m1 m2 m3 }
!
! this is the implementation the formula 8.1.2. (11) from:
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: j1, j2, j3, m1, m2, m3
real(kind=8) :: coeffCG

W3J = 0.0_wp
coeffCG = 0.0_wp
call Clebsh_Gordan(dble(j1)/2.0_wp,dble(m1)/2.0_wp,dble(j2)/2.0_wp,dble(m2)/2.0_wp,dble(j3)/2.0_wp,-dble(m3)/2.0_wp,coeffCG)
if (coeffCG == 0.0_wp) return
W3J = dble((-1)**((j1-j2-m3)/2))*coeffCG/sqrt(dble(j3+1))

return

end function W3J
!=!=
real*8 function WCG(a,al,b,bt,c,gm)
! Calculates a Clebsch-Gordan Coefficient. Argument a, al, b, bt, c, gm are Integer,
! Double their actual value.
! (   c/2, gm/2            )
! ( C                      )
! (   a/2, al/2, b/2, bt/2 )
! this is the implementation the formula 8.2.1. (3) from:
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: a, al, b, bt, c, gm
integer :: lb1, lb2, i
real(kind=8) :: u, fct, dlt
external :: fct, dlt

WCG = 0.0_wp

if ((al+bt) /= gm) return
if (a < 0) return
if (b < 0) return
if (c < 0) return
if (abs(al) > a) return
if (abs(bt) > b) return
if (abs(gm) > c) return
if ((abs(a-b) > c) .or. ((a+b) < c)) return
if ((abs(b-c) > a) .or. ((b+c) < a)) return
if ((abs(c-a) > b) .or. ((c+a) < b)) return
if (mod(a,2) /= mod(abs(al),2)) return
if (mod(b,2) /= mod(abs(bt),2)) return
if (mod(c,2) /= mod(abs(gm),2)) return
u = 0.0_wp
lb1 = min((c-b+al)/2,(c-a-bt)/2)
lb2 = min((a+b-c)/2,(a-al)/2,(b+bt)/2)
if (lb1 < 0) then
  if (-lb1 > lb2) then
    WCG = 0.0_wp
    return
  else
    do i=-lb1,lb2
      u = u+dble((-1)**i)/(fct(i)*fct((a+b-c-2*i)/2)*fct((c-b+al+2*i)/2)*fct((c-a-bt+2*i)/2)*fct((a-al-2*i)/2)*fct((b+bt-2*i)/2))
    end do
  end if
else
  do i=0,lb2
    u = u+dble((-1)**i)/(fct(i)*fct((a+b-c-2*i)/2)*fct((c-b+al+2*i)/2)*fct((c-a-bt+2*i)/2)*fct((a-al-2*i)/2)*fct((b+bt-2*i)/2))
  end do
end if
WCG = u*dlt(a,b,c)*sqrt(fct((a+al)/2)*fct((a-al)/2)*fct((b+bt)/2)*fct((b-bt)/2)*fct((c+gm)/2)*fct((c-gm)/2)*(c+1))

!write(6,'(A)') 'a,  al,  b,  bt,  c,  gm'
!write(6,'(6(F4.1,2x),F20.14)') dble(a)/2.0_wp,dble(al)/2.0_wp,dble(b)/2.0_wp,dble(bt)/2.0_wp,dble(c)/2.0_wp,dble(gm)/2.0_wp,WCG

return

end function WCG
!=!=
real*8 function dlt(a,b,c)
! calculates the delta(a,b,c) function using the formula 8.2.1. from:
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.
!
! a,b,c are positive Integer numbers,
! their values are DoUBLE than their original value

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: a, b, c
real(kind=8) :: fct
logical :: check_triangle
external :: check_triangle, fct

dlt = 0.0_wp
if ((abs(a-b) > c) .or. (a+b < c)) return
if ((abs(b-c) > a) .or. (b+c < a)) return
if ((abs(c-a) > b) .or. (c+a < b)) return
if (mod((a+b-c),2) == 1) return
if (mod((a-b+c),2) == 1) return
if (mod((-a+b+c),2) == 1) return
if (mod((a+b+c),2) == 1) return
if (.not. check_triangle(a,b,c)) return
! special cases:
if (a == 0) dlt = 1.0_wp/sqrt(dble(b+1))
if (b == 0) dlt = 1.0_wp/sqrt(dble(a+1))
if (c == 0) dlt = 1.0_wp/sqrt(dble(a+1))

dlt = sqrt(fct((a+b-c)/2)*fct((a-b+c)/2)*fct((-a+b+c)/2)/fct((a+b+c)/2+1))

return

end function dlt
!=!=
logical function check_triangle(a,b,c)
!  checks If the values a,b,c comply with the triangle rule
implicit none
integer, intent(in) :: a, b, c

check_triangle = .false.

if ((a <= 0) .or. (b <= 0) .or. (c <= 0)) then
  write(6,'(A)') 'a=',a
  write(6,'(A)') 'b=',b
  write(6,'(A)') 'c=',c
  write(6,'(A)') 'The rule is: a>0, b>0 and c>0!'
  write(6,'(A)') 'Please check this issue, or report a bug!'
  return
end if

if ((a+b >= c) .and. (b+c >= a) .and. (c+a >= b)) check_triangle = .true.

return

end function check_triangle
!=!=
complex*16 function WignerD(J,M1,M2,al,bt,gm)
! the function Returns the Wigner-D function specIfying the rotation
! of the |J,M1,M2> around three  angles(alpha,beta,gamma).
! The rotation is either active (ik=1) or passive (ik=2).
!   J, M1, M2 are specIfied as Integer numbers, with a value DoUBLE than their actual size;
!   i.e. J = 2*J (Real); M1= 2*M1(Real); M2=2*M2(Real);
!   alpha, beta, gamma are specIfied as Double precision (Real(kind=8) ::). These values must be defined in
!   radians ( i.e. in units of Pi)
!
!
! This is the implementation of the formula 4.3.(1) from
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: J, M1, M2
real(kind=8), intent(in) :: al, bt, gm
complex(kind=8) :: m1_fact, m2_fact, wig_fac
real(kind=8) :: wigner_d
external :: wigner_d

!  check correctness
WignerD = (0.0_wp,0.0_wp)
wig_fac = (0.0_wp,0.0_wp)
m1_fact = (0.0_wp,0.0_wp)
m2_fact = (0.0_wp,0.0_wp)
if (abs(M1) > J) return
if (abs(M2) > J) return
if (abs(J) < 0) return

m1_fact = exp((0.0_wp,-1.0_wp)**(dble(al*M1)/2.0_wp))
m2_fact = exp((0.0_wp,-1.0_wp)**(dble(gm*M2)/2.0_wp))
wig_fac = cmplx(wigner_d(J,M1,M2,bt),0.0_wp,wp)
WignerD = m1_fact*m2_fact*wig_fac

return

end function WignerD
!=!=
real*8 function wigner_d(J,M1,M2,bt)
! This is the implementation of the formula 4.3.1 (2) from
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: J, M1, M2
real(kind=8), intent(in) :: bt
real(kind=8) :: ksum, fct
integer :: kmin, kmax, i
external :: fct

wigner_d = 0.0_wp
ksum = 0.0_wp
kmax = min((J-M1)/2,(J-M2)/2)
kmin = max(0,-(M1+M2)/2)
do i=kmin,kmax
  ksum = dble((-1)**(i))*(cos(bt/2.0_wp)**dble(M1/2+M2/2+2*i))*(sin(bt/2.0_wp)**dble(J-M1/2-M2/2-2*i))/fct(i)/fct((J-M1)/2-i)/ &
         fct((J-M2)/2-i)/fct((M1+M2)/2+i)
  wigner_d = wigner_d+ksum
end do
wigner_d = wigner_d*dble((-1)**((J-M2)/2))*sqrt(fct((J+M1)/2)*fct((J-M1)/2)*fct((J+M2)/2)*fct((J-M2)/2))

return

end function wigner_d
!=!=
real*8 function RedME(La,Sa,LaP,SaP,L,S)
! function evaluates the reduced matrix elements of the ground atomic J multipet
!
!  the function evaluates the formula:
! <L_A, S_A || Operator ( iL_T, iS_T)  || L_A, S_A >>   formula (S.20) derived by Naoya Iwahara
!
! iL, iS, iL_T, iS_T, Jp, Sp and Lp  are defined as Integers with a value DoUBLE of their true value
!
! the formula is valid for Tb, Dy, Ho, Er, Tm and Yb only

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: La, Sa, LaP, SaP, L, S
integer :: JaP, jm, js, s_orb, l_orb
real(kind=8) :: WCG, temp, factor
external :: WCG

RedME = 0.0_wp
temp = 0.0_wp
l_orb = 6  ! Double of true value l_orb = 3
s_orb = 1  ! Double of true value s_orb = 1/2
JaP = LaP+SaP
if (WCG(La,La,L,0,La,La) == 0.0_wp) return
if (WCG(Sa,Sa,S,0,Sa,Sa) == 0.0_wp) return
if (WCG(La,La,l_orb,LaP-La,LaP,LaP) == 0.0_wp) return
if (WCG(Sa,Sa,s_orb,SaP-Sa,SaP,SaP) == 0.0_wp) return

factor = sqrt(dble((La+1)*(Sa+1)))/WCG(La,La,L,0,La,La)/WCG(Sa,Sa,S,0,Sa,Sa)

do jm=-l_orb,l_orb
  do js=-s_orb,s_orb
    temp = temp+dble((-1)**((l_orb+jm+s_orb+js)/2))*WCG(l_orb,-jm,l_orb,jm,L,0)*WCG(s_orb,-js,s_orb,js,S,0)* &
           WCG(LaP,La+jm,SaP,Sa+js,JaP,La+jm+Sa+js)*WCG(LaP,La+jm,SaP,Sa+js,JaP,La+jm+Sa+js)*WCG(La,La,l_orb,jm,LaP,La+jm)* &
           WCG(La,La,l_orb,jm,LaP,La+jm)*WCG(Sa,Sa,s_orb,js,SaP,Sa+js)*WCG(Sa,Sa,s_orb,js,SaP,Sa+js)/ &
           WCG(La,La,l_orb,LaP-La,LaP,LaP)/WCG(La,La,l_orb,LaP-La,LaP,LaP)/WCG(Sa,Sa,s_orb,SaP-Sa,SaP,SaP)/ &
           WCG(Sa,Sa,s_orb,SaP-Sa,SaP,SaP)
  end do ! is
end do ! im

RedME = temp*factor

return

end function RedME
!=!=
real*8 function jot1(t,L,ML,S,MS,La,Sa,LaP,SaP)
! function evaluates the J1 exchange matrix element ( formula S.19)
!
! iL_T,iL_TM,iS_T,iS_TM,iL,iS,Sp,Lp are defined as Integers with a value DoUBLE of their true value
! the formula is valid for Tb, Dy, Ho, Er, Tm and Yb only
!    Substitutions:

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: L, ML, S, MS, La, Sa, LaP, SaP
integer :: Ja, l_orb
real(kind=8), intent(in) :: t
real(kind=8) :: W9J, WCG, RedME, txt
external :: W9J, WCG, RedME

jot1 = 0.0_wp
txt = 0.0_wp
l_orb = 6  ! Double of true value l_orb = 3
!s_orb = 1 ! Double of true value s_orb = 1/2
if (mod(L,2) == 0) then
  if (ML == 0) then
    txt = dble((-1)**(l_orb/2))*WCG(l_orb,-4,l_orb,4,L,0)
  else if (ML == 8) then
    txt = -0.5_wp*dble((-1)**(l_orb/2))*WCG(l_orb,4,l_orb,4,L,8)
  else if (ML == -8) then
    txt = -0.5_wp*dble((-1)**(l_orb/2))*WCG(l_orb,4,l_orb,4,L,8)
  end if
end if
Ja = La+Sa

jot1 = t*txt*sqrt(dble((Ja+1)*(L+s+1)))*W9J(Ja,La,Sa,Ja,La,Sa,L+s,L,2)*WCG(L,ML,2,MS,L+s,ML+MS)*WCG(Ja,Ja,L+s,0,Ja,Ja)* &
       RedME(La,Sa,LaP,SaP,L,2)/sqrt(2.0_wp)

return

end function jot1
!=!=
real*8 function jot0(t,L,ML,La,Sa,LaP,SaP)
! function evaluates the J1 exchange matrix element ( formula S.19)
!
! iL_T,iL_TM,iS_T,iS_TM,iL,iS,Sp,Lp are defined as Integers with a value DoUBLE of their true value
! the formula is valid for Tb, Dy, Ho, Er, Tm and Yb only
!    Substitutions:

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: L, ML, La, Sa, LaP, SaP
integer :: Ja, l_orb
real(kind=8), intent(in) :: t
real(kind=8) :: W6J, WCG, RedME, txt, W9Jl
external :: W6J, WCG, RedME

jot0 = 0.0_wp
l_orb = 6  ! Double of true value l_orb = 3
!s_orb = 1 ! Double of true value s_orb = 1/2
txt = 0.0_wp
if (mod(L,4) == 0) then
  if (ML == 0) then
    txt = dble((-1)**(l_orb/2))*WCG(l_orb,-4,l_orb,4,L,0)
  else if (ML == 8) then
    txt = -0.5_wp*dble((-1)**(l_orb/2))*WCG(l_orb,4,l_orb,4,L,8)
  else if (ML == -8) then
    txt = -0.5_wp*dble((-1)**(l_orb/2))*WCG(l_orb,4,l_orb,4,L,8)
  end if
end if
Ja = La+Sa
W9Jl = 0.0_wp
W9Jl = dble((-1)**((La+Ja+Sa+L)/2))*sqrt(dble((Sa+1)*(L+1)))*W6J(Ja,La,Sa,La,Ja,L)

jot0 = -t*txt*sqrt(dble((Ja+1)*(L+1)))*W9Jl*WCG(Ja,Ja,L,0,Ja,Ja)*RedME(La,Sa,LaP,SaP,L,0)/sqrt(2.0_wp)

return

end function jot0
!=!=
subroutine verify_CG(N)

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer :: n, k, q, m1, m2
real(wp) :: rJ, rK, rQ, mf, rM1, rM2, rfE, rfG
real(wp) :: CG_A, CG_B, CG_C, CG_D, CG_E, CG_F, CG_G, CG_H
logical :: prn

!2J = n-1
rJ = dble(n-1)/2.0_wp

do k=1,n-1
  do q=0,k
    rK = dble(k)
    rQ = dble(q)

    do m1=1,n
      do m2=1,n
        ! real value for the projections:
        rM1 = -rJ+dble(m1-1)
        rM2 = -rJ+dble(m2-1)

        CG_A = 0.0_wp
        CG_B = 0.0_wp
        CG_C = 0.0_wp
        CG_D = 0.0_wp
        CG_E = 0.0_wp
        CG_F = 0.0_wp
        CG_G = 0.0_wp
        CG_H = 0.0_wp
        mf = (-1)**nint(rK)
        !(a , alpha, b, beta,    c,  gamma)
        call Clebsh_Gordan(rJ,rM2,rK,rQ,rJ,rM1,CG_A)
        call Clebsh_Gordan(rK,rQ,rJ,rM2,rJ,rM1,CG_B)
        call Clebsh_Gordan(rJ,-rM2,rK,-rQ,rJ,-rM1,CG_C)
        call Clebsh_Gordan(rK,-rQ,rJ,-rM2,rJ,-rM1,CG_D)

        rfE = ((-1)**(rJ-rM2))*(sqrt(dble(n)/(2.0_wp*rK+1.0_wp)))
        call Clebsh_Gordan(rJ,rM2,rJ,-rM1,rK,-rQ,CG_E)
        call Clebsh_Gordan(rJ,rM1,rJ,-rM2,rK,rQ,CG_F)

        rfG = (-1)**(rK+rQ)
        call Clebsh_Gordan(rJ,-rM1,rK,rQ,rJ,-rM2,CG_G)
        call Clebsh_Gordan(rK,-rQ,rJ,rM1,rJ,rM2,CG_H)

        prn = (CG_A /= 0.0_wp) .or. (CG_B /= 0.0_wp) .or. (CG_C /= 0.0_wp) .or. (CG_D /= 0.0_wp) .or. (CG_E /= 0.0_wp) .or. &
              (CG_F /= 0.0_wp) .or. (CG_G /= 0.0_wp) .or. (CG_H /= 0.0_wp)

        if (prn) print '(A,1x,8F12.6)','n,k,q,CG:',CG_A,mf*CG_B,mf*CG_C,CG_D,rfE*CG_E,rfE*CG_F,rfG*CG_G,rfG*CG_H

      end do
    end do

  end do
end do

return

end subroutine verify_CG
