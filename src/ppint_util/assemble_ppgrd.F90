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

subroutine Assemble_PPGrd(Fin,nZeta,la,lb,iZeta,Alpha,Beta,A_laplb,A_lamlb,A_lalbp,A_lalbm,JfGrad)

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, iZeta
real(kind=wp), intent(inout) :: Fin(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6)
real(kind=wp), intent(in) :: Alpha, Beta, A_laplb((la+2)*(la+3)/2,(lb+1)*(lb+2)/2), A_lamlb((la+0)*(la+1)/2,(lb+1)*(lb+2)/2), &
                             A_lalbp((la+1)*(la+2)/2,(lb+2)*(lb+3)/2), A_lalbm((la+1)*(la+2)/2,(lb+0)*(lb+1)/2)
logical(kind=iwp), intent(in) :: JfGrad(3,2)
integer(kind=iwp) :: i6, ix, jx, jy, jz

integer(kind=iwp) :: Ind, iy, iz

!                                                                      *
!***********************************************************************
!                                                                      *
! Statement function for cartesian index

Ind(iy,iz) = (iy+iz)*(iy+iz+1)/2+iz+1
!                                                                      *
!***********************************************************************
!                                                                      *
!call RecPrt('X',' ',A_laplb,(la+2)*(la+3)/2,(lb+1)*(lb+2)/2)
!if (la > 0) Call RecPrt('X',' ',A_lamlb,(la+0)*(la+1)/2,(lb+1)*(lb+2)/2)
!call RecPrt('X',' ',A_lalbp,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2)
!if (lb > 0) Call RecPrt('X',' ',A_lalbm,(la+1)*(la+2)/2,(lb+0)*(lb+1)/2)
do ix=la,0,-1
  do iy=la-ix,0,-1
    iz = la-ix-iy

    do jx=lb,0,-1
      do jy=lb-jx,0,-1
        jz = lb-jx-jy

        i6 = 0
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Ax

        if (JfGrad(1,1)) then
          i6 = i6+1
          if (ix == 0) then
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy,iz),Ind(jy,jz))
          else
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy,iz),Ind(jy,jz))- &
                                                  real(ix,kind=wp)*A_lamlb(Ind(iy,iz),Ind(jy,jz))
          end if
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Bx

        if (JfGrad(1,2)) then
          i6 = i6+1
          if (jx == 0) then
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy,jz))
          else
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy,jz))- &
                                                  real(jx,kind=wp)*A_lalbm(Ind(iy,iz),Ind(jy,jz))
          end if
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Ay

        if (JfGrad(2,1)) then
          i6 = i6+1
          if (iy == 0) then
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy+1,iz),Ind(jy,jz))
          else
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy+1,iz),Ind(jy,jz))- &
                                                  real(iy,kind=wp)*A_lamlb(Ind(iy-1,iz),Ind(jy,jz))
          end if
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! By

        if (JfGrad(2,2)) then
          i6 = i6+1
          if (jy == 0) then
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy+1,jz))
          else
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy+1,jz))- &
                                                  real(jy,kind=wp)*A_lalbm(Ind(iy,iz),Ind(jy-1,jz))
          end if
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Az

        if (JfGrad(3,1)) then
          i6 = i6+1
          if (iz == 0) then
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy,iz+1),Ind(jy,jz))
          else
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy,iz+1),Ind(jy,jz))- &
                                                  real(iz,kind=wp)*A_lamlb(Ind(iy,iz-1),Ind(jy,jz))
          end if
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Bz

        if (JfGrad(3,2)) then
          i6 = i6+1
          if (jz == 0) then
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy,jz+1))
          else
            Fin(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy,jz+1))- &
                                                  real(jz,kind=wp)*A_lalbm(Ind(iy,iz),Ind(jy,jz-1))
          end if
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
      end do
    end do
  end do
end do
!call RecPrt('Fin',' ',Fin,nZeta*nElem(la)*nElem(lb),6)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Assemble_PPGrd
