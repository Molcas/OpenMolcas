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

subroutine Assemble_PPGrd(final,nZeta,la,lb,iZeta,Alpha,Beta,A_laplb,A_lamlb,A_lalbp,A_lalbm,JfGrad)

implicit real*8(A-H,O-Z)
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6), A_laplb((la+2)*(la+3)/2,(lb+1)*(lb+2)/2), &
       A_lamlb((la+0)*(la+1)/2,(lb+1)*(lb+2)/2), A_lalbp((la+1)*(la+2)/2,(lb+2)*(lb+3)/2), A_lalbm((la+1)*(la+2)/2,(lb+0)*(lb+1)/2)
logical JfGrad(3,2)

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

        if (.not. JfGrad(1,1)) Go To 102
        i6 = i6+1
        if (ix == 0) then
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy,iz),Ind(jy,jz))
        else
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy,iz),Ind(jy,jz))-dble(ix)*A_lamlb(Ind(iy,iz),Ind(jy,jz))
        end if
102     continue
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Bx

        if (.not. JfGrad(1,2)) Go To 105
        i6 = i6+1
        if (jx == 0) then
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy,jz))
        else
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy,jz))-dble(jx)*A_lalbm(Ind(iy,iz),Ind(jy,jz))
        end if
105     continue
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Ay

        if (.not. JfGrad(2,1)) Go To 103
        i6 = i6+1
        if (iy == 0) then
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy+1,iz),Ind(jy,jz))
        else
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy+1,iz),Ind(jy,jz))- &
                                                  dble(iy)*A_lamlb(Ind(iy-1,iz),Ind(jy,jz))
        end if
103     continue
        !                                                              *
        !***************************************************************
        !                                                              *
        ! By

        if (.not. JfGrad(2,2)) Go To 106
        i6 = i6+1
        if (jy == 0) then
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy+1,jz))
        else
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy+1,jz))- &
                                                  dble(jy)*A_lalbm(Ind(iy,iz),Ind(jy-1,jz))
        end if
106     continue
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Az

        if (.not. JfGrad(3,1)) Go To 104
        i6 = i6+1
        if (iz == 0) then
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy,iz+1),Ind(jy,jz))
        else
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Alpha*A_laplb(Ind(iy,iz+1),Ind(jy,jz))- &
                                                  dble(iz)*A_lamlb(Ind(iy,iz-1),Ind(jy,jz))
        end if
104     continue
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Bz

        if (.not. JfGrad(3,2)) Go To 107
        i6 = i6+1
        if (jz == 0) then
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy,jz+1))
        else
          final(iZeta,Ind(iy,iz),Ind(jy,jz),i6) = Two*Beta*A_lalbp(Ind(iy,iz),Ind(jy,jz+1))- &
                                                  dble(jz)*A_lalbm(Ind(iy,iz),Ind(jy,jz-1))
        end if
107     continue
        !                                                              *
        !***************************************************************
        !                                                              *
      end do
    end do
  end do
end do
!call RecPrt('Final',' ',Final,nZeta*nElem(la)*nElem(lb),6)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Assemble_PPGrd
