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

use wigner_util, only: wcg_real
use Constants, only: Half, cOne
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, k, q
real(kind=wp), intent(out) :: C0
complex(kind=wp), intent(out) :: Cp(n,n), Cm(n,n)
integer(kind=iwp) :: m1, m2
real(kind=wp) :: CGm, CGp, rK, rm1, rm2, rQ, rS
real(kind=wp), external :: fct

rS = real(n-1,kind=wp)*Half
rK = real(k,kind=wp)
rQ = real(q,kind=wp)
C0 = fct(n-1)*sqrt(real(n,kind=wp)/(fct(n-k-1)*fct(n+k)))
! M1 and M2 go from max value to min value
do m1=1,n
  do m2=1,n
    rm1 = rS-real(m1-1,kind=wp)
    rm2 = rS-real(m2-1,kind=wp)
    CGp = wcg_real(rS,rm2,rK,rQ,rS,rm1)
    CGm = wcg_real(rS,rm2,rK,-rQ,rS,rm1)
    Cp(m1,m2) = CGp/C0*cOne
    Cm(m1,m2) = CGm/C0*cOne
  end do
end do

return

end subroutine ITO

subroutine Liviu_ITO(n,k,q,O,W,redME)
! generate O, W as in previous Stewens_matrixel function:
! redME =  ratio between Naoya's and Liviu operators

use Constants, only: cOne
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, k, q
complex(kind=wp), intent(out) :: O(n,n), W(n,n), redME
real(kind=wp) :: C0, CR

call coeff_redus_sub(n,k,CR)
call ITO(n,k,q,C0,O,W)
redME = C0*CR*cOne
W(:,:) = redME*W(:,:)
O(:,:) = redME*O(:,:)

return

end subroutine Liviu_ITO

subroutine Stev_ITO(n,k,q,O,W,redME)
! generate O, W as in previous Stewens_matrixel function:
! redME =  ratio between Naoya's and Liviu operators
! scaled as for Stevens Operators by knm

use Constants, only: One, cZero, cOne
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, k, q
complex(kind=wp), intent(out) :: O(n,n), W(n,n), redME
real(kind=wp) :: F, knm(12,0:12)

O(:,:) = cZero
W(:,:) = cZero
redME = cZero
if ((k > 12) .or. (q > 12)) return

call set_knm(knm)
call Liviu_ITO(n,k,q,O,W,redME)
F = One/knm(k,q)
redME = F*cOne
W(:,:) = redME*W(:,:)
O(:,:) = redME*O(:,:)

return

end subroutine Stev_ITO

subroutine ESO(n,k,q,O,W,redME)
! generate Hermitian ESO operators, as in MATLAB EasySpin(stev) function
! only for q >= 0

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half, cZero, cOne, Onei
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, k, q
complex(kind=wp), intent(out) :: O(n,n), W(n,n), redME
real(kind=wp) :: C0, CR, F, knm(12,0:12)
complex(kind=wp) :: mQ
complex(kind=wp), allocatable :: Cp(:,:), Cm(:,:)

O(:,:) = cZero
W(:,:) = cZero
redME = cZero
if ((k > 12) .or. (q > 12)) return ! not available in MATLAB

call mma_allocate(Cp,n,n,'Cp')
call mma_allocate(Cm,n,n,'Cm')
Cp(:,:) = cZero
Cm(:,:) = cZero

call set_knm(knm)
call coeff_redus_sub(n,k,CR)
call ITO(n,k,q,C0,Cp,Cm)

F = C0*CR/knm(k,q)
redME = F*cOne
mQ = (-cOne)**q
O(:,:) = Half*cOne*redME*(Cm(:,:)+mQ*Cp(:,:))
W(:,:) = Half*Onei*redME*(Cm(:,:)-mQ*Cp(:,:))
call mma_deallocate(Cp)
call mma_deallocate(Cm)

return

end subroutine ESO

subroutine Liviu_ESO(n,k,q,O,W,redME)
! generate Hermitian ESO operators, as in MATLAB EasySpin(stev) function
! only for q >= 0
! Om = Ok-q
! Op = Ok+q
! the only difference with MATLAB's ESO are scaling factor Knm

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half, cOne, Onei
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, k, q
complex(kind=wp), intent(out) :: O(n,n), W(n,n), redME
real(kind=wp) :: C0, CR
complex(kind=wp) :: mQ
complex(kind=wp), allocatable :: Cp(:,:), Cm(:,:)

call mma_allocate(Cp,n,n,'Cp')
call mma_allocate(Cm,n,n,'Cm')

call coeff_redus_sub(n,k,CR)
call ITO(n,k,q,C0,Cp,Cm)
redME = C0*CR*cOne
mQ = (-cOne)**q
O(:,:) = Half*redME*(Cm(:,:)+mQ*Cp(:,:))
W(:,:) = Half*Onei*redME*(Cm(:,:)-mQ*Cp(:,:))
call mma_deallocate(Cp)
call mma_deallocate(Cm)

return

end subroutine Liviu_ESO

subroutine Stewens_matrixel(N,M,d,ITO_O,ITO_W,IPRINT)
! This Subroutine calculates the matrix elements of the ITO  (On)
! on the basis of the eigenfunctions of the effective spin.
!
! N -- the rank of the ITO (On)
! d -- the multiplicity of the effective spin
! Dip_Stewens(N, L, d, d) -- the matrix elements of the ITO tensor
!               operators in the basis of effective spin eigenfunctions

use wigner_util, only: wcg_real
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N, M, d, IPRINT
complex(kind=wp), intent(out) :: ITO_O(d,d), ITO_W(d,d)
integer(kind=iwp) :: i, j, ms1, ms2, npar
real(kind=wp) :: a, al, b, bt, c, COEFF_REDUS, coeffCG, gm
complex(kind=wp), allocatable :: ITO_PLUS(:,:), ITO_MINUS(:,:)
!***********************************************************************

call mma_allocate(ITO_PLUS,[-d,d],[-d,d],label='ITO_PLUS')
call mma_allocate(ITO_MINUS,[-d,d],[-d,d],label='ITO_MINUS')

NPAR = mod(d,2)
ITO_PLUS(:,:) = cZero
ITO_MINUS(:,:) = cZero
ITO_O(:,:) = cZero
ITO_W(:,:) = cZero

call COEFF_REDUS_SUB(d,N,COEFF_REDUS)
!COEFF_REDUS = One
a = real(N,kind=wp)
al = real(M,kind=wp)
c = real(d-1,kind=wp)*Half
b = real(d-1,kind=wp)*Half
do ms1=-(d-NPAR)/2,(d-NPAR)/2
  if ((ms1 == 0) .and. (NPAR == 0)) cycle
  if (NPAR == 0) then
    if (ms1 < 0) then
      gm = real(ms1,kind=wp)+Half
    else
      gm = real(ms1,kind=wp)-Half
    end if
  else
    gm = real(ms1,kind=wp)
  end if

  do ms2=-(d-NPAR)/2,(d-NPAR)/2
    if ((ms2 == 0) .and. (NPAR == 0)) cycle
    if (NPAR == 0) then
      if (ms2 < 0) then
        bt = real(ms2,kind=wp)+Half
      else
        bt = real(ms2,kind=wp)-Half
      end if
    else
      bt = real(ms2,kind=wp)
    end if
    coeffCG = Zero

    coeffCG = wcg_real(a,al,b,bt,c,gm)

    ITO_PLUS(ms1,ms2) = coeffCG*COEFF_REDUS*cOne

    if (iprint > 5) &
      write(u6,'(5x,2(a,i3,2x),6(a,f4.1,2x),2(a,f14.10,2x))') 'ms1=',ms1,'ms2=',ms2,'a=',a,'al=',al,'b=',b,'bt=',bt,'c=',c, &
                                                              'gm=',gm,'coeffCG=',coeffCG,'coeffCG^2=',coeffCG**2
  end do ! ms2
end do ! ms1

al = -real(M,kind=wp)

do ms1=-(d-NPAR)/2,(d-NPAR)/2
  if ((ms1 == 0) .and. (NPAR == 0)) cycle
  if (NPAR == 0) then
    if (ms1 < 0) then
      gm = real(ms1,kind=wp)+Half
    else
      gm = real(ms1,kind=wp)-Half
    end if
  else
    gm = real(ms1,kind=wp)
  end if

  do ms2=-(d-NPAR)/2,(d-NPAR)/2
    if ((ms2 == 0) .and. (NPAR == 0)) cycle

    if (NPAR == 0) then
      if (ms2 < 0) then
        bt = real(ms2,kind=wp)+Half
      else
        bt = real(ms2,kind=wp)-Half
      end if
    else
      bt = real(ms2,kind=wp)
    end if
    coeffCG = Zero

    coeffCG = wcg_real(a,al,b,bt,c,gm)

    ITO_MINUS(ms1,ms2) = coeffCG*COEFF_REDUS*cOne

    if (iprint > 5) &
      write(u6,'(5x,2(a,i3,2x),6(a,f4.1,2x),2(a,f14.10,2x))') 'ms1=',ms1,'ms2=',ms2,'a=',a,'al=',al,'b=',b,'bt=',bt,'c=',c, &
                                                              'gm=',gm,'coeffCG=',coeffCG,'coeffCG^2=',coeffCG**2
  end do
end do

i = 0
do ms1=-(d-NPAR)/2,(d-NPAR)/2
  if ((ms1 == 0) .and. (NPAR == 0)) cycle
  i = i+1
  j = 0

  do ms2=-(d-NPAR)/2,(d-NPAR)/2
    if ((ms2 == 0) .and. (NPAR == 0)) cycle
    j = j+1
    ITO_O(i,j) = ITO_PLUS(ms1,ms2)
    ITO_W(i,j) = ITO_MINUS(ms1,ms2)
  end do
end do

if (IPRINT > 3) then
  write(u6,'(/)')
  write(u6,'(5X,A,2I3)') 'Operator ITO_PLUS',N,M
  write(u6,*)
  do ms1=-(d-NPAR)/2,(d-NPAR)/2
    if ((ms1 == 0) .and. (NPAR == 0)) cycle
    if (NPAR == 1) write(u6,'(16(2X,2ES12.3))') (ITO_PLUS(ms1,ms2),ms2=-(d-NPAR)/2,(d-NPAR)/2)
    if (NPAR == 0) write(u6,'(16(2X,2ES12.3))') (ITO_PLUS(ms1,ms2),ms2=-d/2,-1),(ITO_PLUS(ms1,ms2),ms2=1,d/2)
  end do ! ms1

  write(u6,'(/)')
  write(u6,'(5X,A,2I3)') 'Operator ITO_MINUS',N,M
  write(u6,*)
  do ms1=-(d-NPAR)/2,(d-NPAR)/2
    if ((ms1 == 0) .and. (NPAR == 0)) cycle
    if (NPAR == 1) write(u6,'(16(2X,2ES12.3))') (ITO_MINUS(ms1,ms2),ms2=-(d-NPAR)/2,(d-NPAR)/2)
    if (NPAR == 0) write(u6,'(16(2X,2ES12.3))') (ITO_MINUS(ms1,ms2),ms2=-d/2,-1),(ITO_MINUS(ms1,ms2),ms2=1,d/2)
  end do ! ms1

  write(u6,'(/////)')
  write(u6,'(5X,A,2I3)') 'ITO_O',N,M
  write(u6,*)
  do i=1,d
    write(u6,'(16(2X,2ES12.3))') (ITO_O(i,j),j=1,d)
  end do   ! i
  write(u6,'(/)')
  write(u6,'(5X,A,2I3)') 'ITO_W',N,M
  write(u6,*)
  do i=1,d
    write(u6,'(16(2X,2ES12.3))') (ITO_W(i,j),j=1,d)
  end do ! i
end if ! iPrint

call mma_deallocate(ITO_PLUS)
call mma_deallocate(ITO_MINUS)

return

end subroutine Stewens_matrixel

function fct(n)
! this function provides correct answer till n=169 only

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: fct
integer(kind=iwp), intent(in) :: n
integer(kind=iwp) :: i
real(kind=wp) :: xct, xlim

xct = One
fct = One
if (n < 0) then
  write(u6,'(A,i0)') 'FCT:  N<0 !'
  write(u6,'(A,i0)') 'N = ',N
  write(u6,'(A   )') 'It is an impossible case.'
  fct = -9.0e99_wp
  return

else if (n == 0) then
  return

else if ((n > 0) .and. (n <= 169)) then
  xlim = huge(xct)/real(n,kind=wp)
  do i=1,n
    if (xct > xlim) then
      write(u6,'(A,i0)') 'FCT:   N = ',N
      write(u6,'(A)') 'Factorial of overflows current precision.'
      write(u6,'(A)') 'Use higher numerical precision, or rethink your algorithm.'
      call abend()
    end if
    xct = xct*real(i,kind=wp)
  end do

else
  write(u6,'(A,i0)') 'FCT:   N = ',N
  write(u6,'(A)') 'Factorial of N>169 overflows on x86_64'
  write(u6,'(A)') 'Use higher numerical precision, or rethink your algorithm.'
end if

fct = xct

return

end function fct

subroutine COEFF_REDUS_SUB(d,N,COEFF_REDUS)
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

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: d, N
real(kind=wp), intent(out) :: COEFF_REDUS
real(kind=wp) :: s1, s2
real(kind=wp), parameter :: Norm(51) = [Two**0,Two**1,Two**3,Two**4,Two**7,Two**8,Two**10,Two**11,Two**15,Two**16,Two**18,Two**19, &
                                        Two**22,Two**23,Two**25,Two**26,Two**31,Two**32,Two**34,Two**35,Two**38,Two**39,Two**41, &
                                        Two**42,Two**46,Two**47,Two**49,Two**50,Two**53,Two**54,Two**56,Two**57,Two**63,Two**64, &
                                        Two**66,Two**67,Two**70,Two**71,Two**73,Two**74,Two**78,Zero,Zero,Zero,Zero,Zero,Zero, &
                                        Zero,Zero,Zero,Zero]
real(kind=wp), external :: fct

!  new method
s1 = sqrt(fct(d+N)/fct(d-N-1))
s2 = real(2**N,kind=wp)*sqrt(real(d,kind=wp))

COEFF_REDUS = Norm(N/2+1)*s1/s2

return

end subroutine COEFF_REDUS_SUB

function dlt(a,b,c)
! calculates the delta(a,b,c) function using the formula 8.2.1. from:
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World Scientific, 1988.
!
! a,b,c are positive Integer numbers,
! their values are DoUBLE than their original value

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: dlt
integer(kind=iwp), intent(in) :: a, b, c
real(kind=wp), external :: fct
logical(kind=iwp), external :: check_triangle

dlt = Zero
if ((abs(a-b) > c) .or. (a+b < c)) return
if ((abs(b-c) > a) .or. (b+c < a)) return
if ((abs(c-a) > b) .or. (c+a < b)) return
if (mod((a+b-c),2) == 1) return
if (mod((a-b+c),2) == 1) return
if (mod((-a+b+c),2) == 1) return
if (mod((a+b+c),2) == 1) return
if (.not. check_triangle(a,b,c)) return
! special cases:
if (a == 0) dlt = One/sqrt(real(b+1,kind=wp))
if (b == 0) dlt = One/sqrt(real(a+1,kind=wp))
if (c == 0) dlt = One/sqrt(real(a+1,kind=wp))

dlt = sqrt(fct((a+b-c)/2)*fct((a-b+c)/2)*fct((-a+b+c)/2)/fct((a+b+c)/2+1))

return

end function dlt

function check_triangle(a,b,c)
!  checks If the values a,b,c comply with the triangle rule

use Definitions, only: iwp, u6

implicit none
logical(kind=iwp) :: check_triangle
integer(kind=iwp), intent(in) :: a, b, c

check_triangle = .false.

if ((a <= 0) .or. (b <= 0) .or. (c <= 0)) then
  write(u6,'(A)') 'a=',a
  write(u6,'(A)') 'b=',b
  write(u6,'(A)') 'c=',c
  write(u6,'(A)') 'The rule is: a>0, b>0 and c>0!'
  write(u6,'(A)') 'Please check this issue, or report a bug!'
  return
end if

if ((a+b >= c) .and. (b+c >= a) .and. (c+a >= b)) check_triangle = .true.

return

end function check_triangle

function WignerD(J,M1,M2,al,bt,gm)
! the function Returns the Wigner-D function specifying the rotation
! of the |J,M1,M2> around three  angles(alpha,beta,gamma).
! The rotation is either active (ik=1) or passive (ik=2).
!   J, M1, M2 are specified as Integer numbers, with a value DoUBLE than their actual size;
!   i.e. J = 2*J (Real); M1= 2*M1(Real); M2=2*M2(Real);
!   alpha, beta, gamma are specified as double precision (Real(kind=wp) ::). These values must be defined in
!   radians ( i.e. in units of Pi)
!
!
! This is the implementation of the formula 4.3.(1) from
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World Scientific, 1988.

use Constants, only: Half, cZero, cOne, Onei
use Definitions, only: wp, iwp

implicit none
complex(kind=wp) :: WignerD
integer(kind=iwp), intent(in) :: J, M1, M2
real(kind=wp), intent(in) :: al, bt, gm
complex(kind=wp) :: m1_fact, m2_fact, wig_fac
real(kind=wp), external :: wigner_d

! check correctness
WignerD = cZero
if (abs(M1) > J) return
if (abs(M2) > J) return
if (abs(J) < 0) return

m1_fact = exp((-Onei)**(real(al*M1,kind=wp)*Half))
m2_fact = exp((-Onei)**(real(gm*M2,kind=wp)*Half))
wig_fac = wigner_d(J,M1,M2,bt)*cOne
WignerD = m1_fact*m2_fact*wig_fac

return

end function WignerD

function wigner_d(J,M1,M2,bt)
! This is the implementation of the formula 4.3.1 (2) from
!   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
!   "Quantum Theory of Angular Momentum", World Scientific, 1988.

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: wigner_d
integer(kind=iwp), intent(in) :: J, M1, M2
real(kind=wp), intent(in) :: bt
integer(kind=iwp) :: i, kmax, kmin
real(kind=wp) :: ksum
real(kind=wp), external :: fct

wigner_d = Zero
ksum = Zero
kmax = min((J-M1)/2,(J-M2)/2)
kmin = max(0,-(M1+M2)/2)
do i=kmin,kmax
  ksum = real((-1)**i,kind=wp)*(cos(bt*Half)**(M1/2+M2/2+2*i))*(sin(bt*Half)**(J-M1/2-M2/2-2*i))/fct(i)/fct((J-M1)/2-i)/ &
         fct((J-M2)/2-i)/fct((M1+M2)/2+i)
  wigner_d = wigner_d+ksum
end do
wigner_d = wigner_d*real((-1)**((J-M2)/2),kind=wp)*sqrt(fct((J+M1)/2)*fct((J-M1)/2)*fct((J+M2)/2)*fct((J-M2)/2))

return

end function wigner_d

function RedME(La,Sa,LaP,SaP,L,S)
! function evaluates the reduced matrix elements of the ground atomic J multipet
!
!  the function evaluates the formula:
! <L_A, S_A || Operator ( iL_T, iS_T)  || L_A, S_A >>   formula (S.20) derived by Naoya Iwahara
!
! iL, iS, iL_T, iS_T, Jp, Sp and Lp  are defined as Integers with a value DoUBLE of their true value
!
! the formula is valid for Tb, Dy, Ho, Er, Tm and Yb only

use wigner_util, only: wcg
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: RedME
integer(kind=iwp), intent(in) :: La, Sa, LaP, SaP, L, S
integer(kind=iwp) :: JaP, jm, js, l_orb, s_orb
real(kind=wp) :: factor, temp

RedME = Zero
l_orb = 6  ! double of true value l_orb = 3
s_orb = 1  ! double of true value s_orb = 1/2
JaP = LaP+SaP
if (WCG(La,La,L,0,La,La) == Zero) return
if (WCG(Sa,Sa,S,0,Sa,Sa) == Zero) return
if (WCG(La,La,l_orb,LaP-La,LaP,LaP) == Zero) return
if (WCG(Sa,Sa,s_orb,SaP-Sa,SaP,SaP) == Zero) return

factor = sqrt(real((La+1)*(Sa+1),kind=wp))/WCG(La,La,L,0,La,La)/WCG(Sa,Sa,S,0,Sa,Sa)

temp = Zero
do jm=-l_orb,l_orb
  do js=-s_orb,s_orb
    temp = temp+real((-1)**((l_orb+jm+s_orb+js)/2),kind=wp)*WCG(l_orb,-jm,l_orb,jm,L,0)*WCG(s_orb,-js,s_orb,js,S,0)* &
           WCG(LaP,La+jm,SaP,Sa+js,JaP,La+jm+Sa+js)*WCG(LaP,La+jm,SaP,Sa+js,JaP,La+jm+Sa+js)*WCG(La,La,l_orb,jm,LaP,La+jm)* &
           WCG(La,La,l_orb,jm,LaP,La+jm)*WCG(Sa,Sa,s_orb,js,SaP,Sa+js)*WCG(Sa,Sa,s_orb,js,SaP,Sa+js)/ &
           WCG(La,La,l_orb,LaP-La,LaP,LaP)/WCG(La,La,l_orb,LaP-La,LaP,LaP)/WCG(Sa,Sa,s_orb,SaP-Sa,SaP,SaP)/ &
           WCG(Sa,Sa,s_orb,SaP-Sa,SaP,SaP)
  end do ! is
end do ! im

RedME = temp*factor

return

end function RedME

function jot1(t,L,ML,S,MS,La,Sa,LaP,SaP)
! function evaluates the J1 exchange matrix element ( formula S.19)
!
! iL_T,iL_TM,iS_T,iS_TM,iL,iS,Sp,Lp are defined as Integers with a value DoUBLE of their true value
! the formula is valid for Tb, Dy, Ho, Er, Tm and Yb only
!    Substitutions:

use wigner_util, only: w9j, wcg
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: jot1
real(kind=wp), intent(in) :: t
integer(kind=iwp), intent(in) :: L, ML, S, MS, La, Sa, LaP, SaP
integer(kind=iwp) :: Ja, l_orb
real(kind=wp) :: txt
real(kind=wp), external :: RedME

jot1 = Zero
l_orb = 6  ! double of true value l_orb = 3
!s_orb = 1 ! double of true value s_orb = 1/2
txt = Zero
if (mod(L,2) == 0) then
  if (ML == 0) then
    txt = real((-1)**(l_orb/2),kind=wp)*WCG(l_orb,-4,l_orb,4,L,0)
  else if (ML == 8) then
    txt = -Half*real((-1)**(l_orb/2),kind=wp)*WCG(l_orb,4,l_orb,4,L,8)
  else if (ML == -8) then
    txt = -Half*real((-1)**(l_orb/2),kind=wp)*WCG(l_orb,4,l_orb,4,L,8)
  end if
end if
Ja = La+Sa

jot1 = t*txt*sqrt(real((Ja+1)*(L+s+1),kind=wp))*W9J(Ja,La,Sa,Ja,La,Sa,L+s,L,2)*WCG(L,ML,2,MS,L+s,ML+MS)*WCG(Ja,Ja,L+s,0,Ja,Ja)* &
       RedME(La,Sa,LaP,SaP,L,2)/sqrt(Two)

return

end function jot1

function jot0(t,L,ML,La,Sa,LaP,SaP)
! function evaluates the J1 exchange matrix element ( formula S.19)
!
! iL_T,iL_TM,iS_T,iS_TM,iL,iS,Sp,Lp are defined as Integers with a value DoUBLE of their true value
! the formula is valid for Tb, Dy, Ho, Er, Tm and Yb only
!    Substitutions:

use wigner_util, only: w6j, wcg
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: jot0
real(kind=wp), intent(in) :: t
integer(kind=iwp), intent(in) :: L, ML, La, Sa, LaP, SaP
integer(kind=iwp) :: Ja, l_orb
real(kind=wp) :: txt, W9Jl
real(kind=wp), external :: RedME

jot0 = Zero
l_orb = 6  ! double of true value l_orb = 3
!s_orb = 1 ! double of true value s_orb = 1/2
txt = Zero
if (mod(L,4) == 0) then
  if (ML == 0) then
    txt = real((-1)**(l_orb/2),kind=wp)*WCG(l_orb,-4,l_orb,4,L,0)
  else if (ML == 8) then
    txt = -Half*real((-1)**(l_orb/2),kind=wp)*WCG(l_orb,4,l_orb,4,L,8)
  else if (ML == -8) then
    txt = -Half*real((-1)**(l_orb/2),kind=wp)*WCG(l_orb,4,l_orb,4,L,8)
  end if
end if
Ja = La+Sa
W9Jl = real((-1)**((La+Ja+Sa+L)/2),kind=wp)*sqrt(real((Sa+1)*(L+1),kind=wp))*W6J(Ja,La,Sa,La,Ja,L)

jot0 = -t*txt*sqrt(real((Ja+1)*(L+1),kind=wp))*W9Jl*WCG(Ja,Ja,L,0,Ja,Ja)*RedME(La,Sa,LaP,SaP,L,0)/sqrt(Two)

return

end function jot0

subroutine verify_CG(N)

use wigner_util, only: wcg_real
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
integer(kind=iwp) :: k, m1, m2, q
real(kind=wp) :: CG_A, CG_B, CG_C, CG_D, CG_E, CG_F, CG_G, CG_H, mf, rfE, rfG, rJ, rK, rM1, rM2, rQ
logical(kind=iwp) :: prn

!2J = n-1
rJ = real(n-1,kind=wp)*Half

do k=1,n-1
  do q=0,k
    rK = real(k,kind=wp)
    rQ = real(q,kind=wp)

    do m1=1,n
      do m2=1,n
        ! real value for the projections:
        rM1 = -rJ+real(m1-1,kind=wp)
        rM2 = -rJ+real(m2-1,kind=wp)

        mf = (-1)**nint(rK)
        ! (a , alpha, b, beta,    c,  gamma)
        CG_A = wcg_real(rJ,rM2,rK,rQ,rJ,rM1)
        CG_B = wcg_real(rK,rQ,rJ,rM2,rJ,rM1)
        CG_C = wcg_real(rJ,-rM2,rK,-rQ,rJ,-rM1)
        CG_D = wcg_real(rK,-rQ,rJ,-rM2,rJ,-rM1)

        rfE = ((-1)**(rJ-rM2))*(sqrt(real(n,kind=wp)/(Two*rK+One)))
        CG_E = wcg_real(rJ,rM2,rJ,-rM1,rK,-rQ)
        CG_F = wcg_real(rJ,rM1,rJ,-rM2,rK,rQ)

        rfG = (-1)**(rK+rQ)
        CG_G = wcg_real(rJ,-rM1,rK,rQ,rJ,-rM2)
        CG_H = wcg_real(rK,-rQ,rJ,rM1,rJ,rM2)

        prn = (CG_A /= Zero) .or. (CG_B /= Zero) .or. (CG_C /= Zero) .or. (CG_D /= Zero) .or. (CG_E /= Zero) .or. &
              (CG_F /= Zero) .or. (CG_G /= Zero) .or. (CG_H /= Zero)

        if (prn) print '(A,1x,8F12.6)','n,k,q,CG:',CG_A,mf*CG_B,mf*CG_C,CG_D,rfE*CG_E,rfE*CG_F,rfG*CG_G,rfG*CG_H

      end do
    end do

  end do
end do

return

end subroutine verify_CG
