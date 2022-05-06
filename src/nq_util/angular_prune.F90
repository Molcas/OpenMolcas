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

subroutine Angular_Prune(Radius,nR,iAngular_Grid,Crowding,Fade,R_BS,L_Quad,R_Min,lAng,nAngularGrids)

use nq_Structure, only: Info_Ang
implicit real*8(a-h,o-z)
real*8 Radius(2,nR), R_Min(0:lAng)
integer iAngular_Grid(nR)
#include "real.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
R_Test = R_BS/Crowding
#ifdef _DEBUGPRINT_
write(6,*) 'lAng=',lAng
write(6,*) 'Crowding=',Crowding
write(6,*) 'R_BS=',R_BS
write(6,*) 'L_Quad=',L_Quad
write(6,*) 'LMax_NQ=',size(Info_Ang)
write(6,*) 'nAngularGrids=',nAngularGrids
write(6,*) 'R_Test=',R_Test
write(6,'(A,10G10.3)') 'R_Min=',R_Min
write(6,*) 'Info_Ang(*)%L_Eff=',(Info_Ang(i)%L_Eff,i=1,nAngularGrids)
#endif
do iR=1,nR

  R_Value = Radius(1,iR)
# ifdef _DEBUGPRINT_
  write(6,'(A,G10.3)') 'R_Value=',R_Value
# endif

  ! Establish L_Eff according to the crowding factor.

  ! Avoid overflow by converting to Int at the end
  iAng = int(Half*min(L_Quad*R_Value/R_Test,dble(L_Quad)))
# ifdef _DEBUGPRINT_
  write(6,*) 'iAng=',iAng
# endif

  ! Close to the nuclei we can use the alternative formula
  ! given by the LMG radial grid.

  iAng = max(iAng,lAng)
  do jAng=lAng,1,-1
    if (R_Value < R_Min(jAng)) iAng = min(iAng,jAng-1)
  end do
# ifdef _DEBUGPRINT_
  write(6,*) 'iAng=',iAng
# endif

  ! Fade the outer part

  R_Test2 = Fade*R_BS
  ! Avoid overflow by converting to Int at the end
  if (R_Value > R_Test2) iAng = int(Half*min(L_Quad*R_Test2/R_Value,dble(L_Quad)))

  ! Since the Lebedev grid is not defined for any L_Eff
  ! value we have to find the closest one, i.e. of the same
  ! order or higher. Start loop from low order!

  kSet = 0
  do jSet=1,nAngularGrids
    if ((Info_Ang(jSet)%L_Eff >= 2*iAng+1) .and. (kSet == 0)) then
      kSet = jSet
      iAng = Info_Ang(kSet)%L_Eff/2
    end if
  end do
  if (kSet == 0) then
    kSet = nAngularGrids
    iAng = Info_Ang(kSet)%L_Eff/2
  end if

  iAngular_Grid(iR) = kSet

end do
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'iAngular_Grid:'
write(6,*) 'R_Min     R_Max    kSet  nR   nPoints'
RStart = Radius(1,1)
kSet = iAngular_Grid(1)
nTot = Info_Ang(1)%nPoints
mR = 1
do iR=2,nR
  nTot = nTot+Info_Ang(iAngular_Grid(iR))%nPoints
  if ((iAngular_Grid(iR) /= kSet) .or. (iR == nR)) then
    jR = iR-1
    if (iR == nR) then
      jR = nR
      mR = mR+1
    end if
    write(6,'(2G10.3,I3,2I5)') RStart,Radius(1,jR),kSet,mR,Info_Ang(kSet)%nPoints
    RStart = Radius(1,iR)
    kSet = iAngular_Grid(iR)
    mR = 0
  end if
  mR = mR+1
end do
write(6,*) 'Total grid size:',nTot
write(6,*)
#endif

return

end subroutine Angular_Prune
