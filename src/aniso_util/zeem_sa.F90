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

use Index_Functions, only: iTri, nTri_Elem
use Constants, only: Zero, cZero, cOne, cm_s, hPlanck, mBohr
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: H, dX, dY, dZ, W(N), S(3), zJ
complex(kind=wp), intent(in) :: M(3,N,N), sM(3,N,N)
real(kind=wp), intent(out) :: WM(N), RWORK(3*N-2)
complex(kind=wp), intent(out) :: ZM(N,N), HZEE(nTri_Elem(N)), WORK(2*N-1), W_c(N)
logical(kind=iwp), intent(in) :: DBG
integer(kind=iwp) :: i, info, j
complex(kind=wp) :: dX_c, dY_c, dZ_c, H_c, mB_c, P, R, RP, S_c(3), zJ_c
real(kind=wp), parameter :: mB = mBohr/(cm_s*hPlanck) ! in cm-1*T-1

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

      HZEE(iTri(i,j)) = HZEE(iTri(i,j))-mB_c*H_c*R

    end do
  end do

else ! zJ /= 0

  do i=1,N
    do j=1,i
      R = dX_c*M(1,j,i)+dY_c*M(2,j,i)+dZ_c*M(3,j,i)

      P = dX_c*SM(1,j,i)*S_c(1)+dY_c*SM(2,j,i)*S_c(2)+dZ_c*SM(3,j,i)*S_c(3)

      RP = mB_c*H_c*R+zJ_c*P

      HZEE(iTri(i,j)) = HZEE(iTri(i,j))-RP
    end do
  end do

end if ! zJ

! add diagonal energies:
do i=1,N
  HZEE(iTri(i,i)) = HZEE(iTri(i,i))+W_c(i)
end do

if (DBG) then
  write(u6,'(A)') 'HZEE:'
  do i=1,N
    do j=1,i
      write(u6,'(2i3,2x,100(2ES16.8,2x))') i,j,HZEE(iTri(i,j))
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
