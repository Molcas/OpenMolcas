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

implicit none
integer, parameter :: wp = kind(0.d0)
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
real(kind=8) :: mB
real(kind=8) :: RWORK(3*N-2)
complex(kind=8) :: HZEE(N*(N+1)/2)
complex(kind=8) :: WORK(2*N-1), R, P, RP
complex(kind=8) :: H_c, dX_c, dY_c, dZ_c, zJ_c, W_c(N), S_c(3)
complex(kind=8) :: mB_c
logical :: DBG

mB = 0.4668643740_wp ! in cm-1*T-1

! initialization

!call dcopy_(N,[0.0_wp],0,WM,1)
!call dcopy_(3*N-2,[0.0_wp],0,RWORK,1)
!call zcopy_(N**2,[(0.0_wp,0.0_wp)],0,ZM,1)
!call zcopy_(N*(N+1)/2,[(0.0_wp,0.0_wp)],0,HZEE,1)
!call zcopy_(2*N-1,[(0.0_wp,0.0_wp)],0,WORK,1)
WM = 0.0_wp
RWORK = 0.0_wp
ZM = (0.0_wp,0.0_wp)
HZEE = (0.0_wp,0.0_wp)
WORK = (0.0_wp,0.0_wp)

info = 0
if (DBG) then
  write(6,'(A,4ES20.10)') 'dX,dY,dZ,H =',dX,dY,dZ,H
  write(6,'(A,4ES20.10)') 'Sx,y,z=',S(1),S(2),S(3)
  write(6,*) 'zJ = ',zJ
end if

! initialize
H_c = cmplx(H,0.0_wp,wp)
dX_c = cmplx(dX,0.0_wp,wp)
dY_c = cmplx(dY,0.0_wp,wp)
dZ_c = cmplx(dZ,0.0_wp,wp)
zJ_c = cmplx(zJ,0.0_wp,wp)
mB_c = cmplx(mB,0.0_wp,wp)
call zcopy_(N,[(0.0_wp,0.0_wp)],0,W_c,1)
call zcopy_(3,[(0.0_wp,0.0_wp)],0,S_c,1)
do i=1,N
  W_c(i) = cmplx(W(i),0.0_wp,wp)
end do
do i=1,3
  S_c(i) = cmplx(S(i),0.0_wp,wp)
end do

if (DBG) then
  write(6,*) ' H_c = ',H_c
  write(6,*) 'dX_c = ',dX_c
  write(6,*) 'dY_c = ',dY_c
  write(6,*) 'dZ_c = ',dZ_c
  write(6,*) 'zJ_c = ',zJ_c
  write(6,*) 'mB_c = ',mB_c
end if

! build the Zeeman Hamiltonian
if (abs(zJ) < tiny(0.0_wp)) then
  ! zJ = 0

  do i=1,N
    do j=1,i
      R = (0.0_wp,0.0_wp)
      R = dX_c*M(1,j,i)+dY_c*M(2,j,i)+dZ_c*M(3,j,i)

      HZEE(j+(i-1)*i/2) = HZEE(j+(i-1)*i/2)-mB_c*H_c*R

    end do
  end do

else ! zJ /= 0

  do i=1,N
    do j=1,i
      R = (0.0_wp,0.0_wp)
      P = (0.0_wp,0.0_wp)
      RP = (0.0_wp,0.0_wp)

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
  write(6,'(A)') 'HZEE:'
  do i=1,N
    do j=1,i
      write(6,'(2i3,2x,100(2ES16.8,2x))') i,j,HZEE(j+(i-1)*i/2)
    end do
  end do
end if
! diagonalization
call zhpev_('V','U',N,HZEE,WM,ZM,N,WORK,RWORK,INFO)

if (DBG) then
  do i=1,N
    write(6,'(A,i3,A,F20.13,A,i2,A,99(2F16.10,1x))') 'WM(',i,')=',WM(i)!,' ZM(j,',i,'):',(ZM(j,i),j=1,N)
  end do
end if

return

end subroutine ZEEM_SA
