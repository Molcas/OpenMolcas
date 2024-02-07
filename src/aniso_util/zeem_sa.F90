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

subroutine ZEEM_SA(N,H,dX,dY,dZ,W,M,sM,S,zJ,WM,ZM,DBG,RWORK,HZEE,WORK,W_c)

use Constants, only: Zero, cZero, cOne, cLight, mBohr, rPlanck
use Definitions, only: wp, u6

implicit none
! input variables:
integer, intent(in) :: N
real(kind=8), intent(in) :: H, dX, dY, dZ, zJ
real(kind=8), intent(in) :: W(N)
real(kind=8), intent(in) :: S(3)
complex(kind=8), intent(in) :: sM(3,N,N)
complex(kind=8), intent(in) :: M(3,N,N)
! output variables:
real(kind=8), intent(out) :: WM(N)
complex(kind=8), intent(out) :: ZM(N,N)
! local variables:
integer :: i, j, info
real(kind=8), parameter :: mB = mBohr/(cLight*rPlanck*1.0e2_wp) ! in cm-1*T-1
real(kind=8) :: RWORK(3*N-2)
complex(kind=8) :: HZEE(N*(N+1)/2)
complex(kind=8) :: WORK(2*N-1), R, P, RP
complex(kind=8) :: H_c, dX_c, dY_c, dZ_c, zJ_c, W_c(N), S_c(3)
complex(kind=8) :: mB_c
logical :: DBG

! initialization

WM(:) = Zero
RWORK(:) = Zero
ZM(:,:) = cZero
HZEE(:) = cZero
WORK(:) = cZero

info = 0
if (DBG) then
  write(u6,'(A,4ES20.10)') 'dX,dY,dZ,H =',dX,dY,dZ,H
  write(u6,'(A,4ES20.10)') 'Sx,y,z=',S(1),S(2),S(3)
  write(u6,*) 'zJ = ',zJ
end if

! initialize
H_c = H*cOne
dX_c = dX*cOne
dY_c = dY*cOne
dZ_c = dZ*cOne
zJ_c = zJ*cOne
mB_c = mB*cOne
W_c(:) = W(:)*cOne
S_c(:) = S(:)*cOne

if (DBG) then
  write(u6,*) ' H_c = ',H_c
  write(u6,*) 'dX_c = ',dX_c
  write(u6,*) 'dY_c = ',dY_c
  write(u6,*) 'dZ_c = ',dZ_c
  write(u6,*) 'zJ_c = ',zJ_c
  write(u6,*) 'mB_c = ',mB_c
end if

! build the Zeeman Hamiltonian
if (abs(zJ) < tiny(zJ)) then
  ! zJ = 0

  do i=1,N
    do j=1,i
      R = dX_c*M(1,j,i)+dY_c*M(2,j,i)+dZ_c*M(3,j,i)

      HZEE(j+(i-1)*i/2) = HZEE(j+(i-1)*i/2)-mB_c*H_c*R

    end do
  end do

else ! zJ /= 0

  do i=1,N
    do j=1,i
      R = dX_c*M(1,j,i)+dY_c*M(2,j,i)+dZ_c*M(3,j,i)

      P = dX_c*SM(1,j,i)*S_c(1)+dY_c*SM(2,j,i)*S_c(2)+dZ_c*SM(3,j,i)*S_c(3)

      RP = mB_c*H_c*R+zJ_c*P

      HZEE(j+(i-1)*i/2) = HZEE(j+(i-1)*i/2)-RP
    end do
  end do

end if ! zJ

! add diagonal energies:
do i=1,N
  HZEE(i+(i-1)*i/2) = HZEE(i+(i-1)*i/2)+W_c(i)
end do

if (DBG) then
  write(u6,'(A)') 'HZEE:'
  do i=1,N
    do j=1,i
      write(u6,'(2i3,2x,100(2ES16.8,2x))') i,j,HZEE(j+(i-1)*i/2)
    end do
  end do
end if
! diagonalization
call zhpev_('V','U',N,HZEE,WM,ZM,N,WORK,RWORK,INFO)

if (DBG) then
  do i=1,N
    write(u6,'(A,i3,A,F20.13,A,i2,A,99(2F16.10,1x))') 'WM(',i,')=',WM(i)!,' ZM(j,',i,'):',(ZM(j,i),j=1,N)
  end do
end if

return

end subroutine ZEEM_SA
