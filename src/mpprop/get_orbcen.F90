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

subroutine Get_OrbCen(nPrim,NORBI,Q_MltPl,RCHC,CENTX,CENTY,CENTZ,OCOF)
!EB subroutine Get_OrbCen(nPrim,nBas,NORBI,Q_MltPl,RCHC,

implicit real*8(a-h,o-z)

dimension RCPO(3,NORBI), RCMI(3,NORBI)
dimension CENTX(nPrim*(nPrim+1)/2)
dimension CENTY(nPrim*(nPrim+1)/2)
dimension CENTZ(nPrim*(nPrim+1)/2)
dimension Q_MltPl(nPrim*(nPrim+1)/2)
dimension RCHC(3,NORBI)
dimension oCof(NORBI,nPrim)

! CALCULATE A CENTER OF CHARGE FOR EACH MOLECULAR ORBITAL

do I=1,NORBI
  RCPO(1,I) = 0.0
  RCPO(2,I) = 0.0
  RCPO(3,I) = 0.0
  RCMI(1,I) = 0.0
  RCMI(2,I) = 0.0
  RCMI(3,I) = 0.0
  QPOS = 0.0
  QMIN = 0.0
  do J=1,nPrim
    do K=1,J

      ! THE WEIGHTED AVERAGE OF THE POSITIVE AND NEGATIVE CONTRIBUTIONS

      OOQ = OCOF(I,J)*OCOF(I,K)*Q_MltPl(J*(J-1)/2+K)*2.0
      if (OOQ >= 0.0) then
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
    if (OOQ >= 0.0) then
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
  RCHC(1,I) = (RCPO(1,I)-RCMI(1,I))/(QPOS-QMIN)
  RCHC(2,I) = (RCPO(2,I)-RCMI(2,I))/(QPOS-QMIN)
  RCHC(3,I) = (RCPO(3,I)-RCMI(3,I))/(QPOS-QMIN)
end do

return

!EB 96 format(I5,3E15.8)

end subroutine Get_OrbCen
