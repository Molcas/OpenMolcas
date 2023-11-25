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

subroutine Get_OrbCen(nPrim,NORBI,Q_MltPl,RCHC,CENTX,CENTY,CENTZ,oCof)
!EB subroutine Get_OrbCen(nPrim,nBas,NORBI,Q_MltPl,RCHC,

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nPrim, NORBI
real(kind=wp), intent(in) :: Q_MltPl(nPrim*(nPrim+1)/2), CENTX(nPrim*(nPrim+1)/2), CENTY(nPrim*(nPrim+1)/2), &
                             CENTZ(nPrim*(nPrim+1)/2), oCof(NORBI,nPrim)
real(kind=wp), intent(out) :: RCHC(3,NORBI)
integer(kind=iwp) :: I, J, K
real(kind=wp) :: OOQ, QMIN, QPOS
real(kind=wp), allocatable :: RCPO(:,:), RCMI(:,:)

! CALCULATE A CENTER OF CHARGE FOR EACH MOLECULAR ORBITAL

call mma_allocate(RCPO,3,NORBI,label='RCPO')
call mma_allocate(RCMI,3,NORBI,label='RCMI')

RCPO(:,:) = Zero
RCMI(:,:) = Zero
do I=1,NORBI
  QPOS = Zero
  QMIN = Zero
  do J=1,nPrim
    do K=1,J

      ! THE WEIGHTED AVERAGE OF THE POSITIVE AND NEGATIVE CONTRIBUTIONS

      OOQ = OCOF(I,J)*OCOF(I,K)*Q_MltPl(J*(J-1)/2+K)*Two
      if (OOQ >= Zero) then
        QPOS = QPOS+OOQ
        RCPO(1,I) = RCPO(1,I)+OOQ*CENTX(J*(J-1)/2+K)
        RCPO(2,I) = RCPO(2,I)+OOQ*CENTY(J*(J-1)/2+K)
        RCPO(3,I) = RCPO(3,I)+OOQ*CENTZ(J*(J-1)/2+K)
      else
        QMIN = QMIN+OOQ
        RCMI(1,I) = RCMI(1,I)+OOQ*CENTX(J*(J-1)/2+K)
        RCMI(2,I) = RCMI(2,I)+OOQ*CENTY(J*(J-1)/2+K)
        RCMI(3,I) = RCMI(3,I)+OOQ*CENTZ(J*(J-1)/2+K)
      end if
    end do
    OOQ = OCOF(I,J)*OCOF(I,J)*Q_MltPl(J*(J+1)/2)
    if (OOQ >= Zero) then
      QPOS = QPOS-OOQ
      RCPO(1,I) = RCPO(1,I)-OOQ*CENTX(J*(J+1)/2)
      RCPO(2,I) = RCPO(2,I)-OOQ*CENTY(J*(J+1)/2)
      RCPO(3,I) = RCPO(3,I)-OOQ*CENTZ(J*(J+1)/2)
    else
      QMIN = QMIN-OOQ
      RCMI(1,I) = RCMI(1,I)-OOQ*CENTX(J*(J+1)/2)
      RCMI(2,I) = RCMI(2,I)-OOQ*CENTY(J*(J+1)/2)
      RCMI(3,I) = RCMI(3,I)-OOQ*CENTZ(J*(J+1)/2)
    end if
  end do
  RCHC(:,I) = (RCPO(:,I)-RCMI(:,I))/(QPOS-QMIN)
end do

call mma_deallocate(RCPO)
call mma_deallocate(RCMI)

return

!EB 96 format(I5,3ES15.8)

end subroutine Get_OrbCen
