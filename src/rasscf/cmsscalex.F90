!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 12, 2022, created this file.               *
! ****************************************************************

subroutine CMSScaleX(X,R,DeltaR,Qnew,Qold,RCopy,GDCopy,DgCopy,GDstate,GDOrbit,Dgstate,DgOrbit,DDg,nSPair,lRoots2,nGD,NAC2,nDDg, &
                     Saved)

use CMS, only: NCMSScale
use rasscf_global, only: CMSThreshold, lRoots
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nSPair, lRoots2, nGD, NAC2, nDDg
real(kind=wp) :: X(nSPair), R(lRoots2), DeltaR(lRoots2), Qnew, QOld, RCopy(lRoots2), GDCopy(nGD), DgCopy(nGD), GDState(nGD), &
                 GDOrbit(nGD), Dgstate(nGD), DgOrbit(nGD), DDg(nDDg)
logical(kind=iwp) :: Saved
integer(kind=iwp) :: nScaleMax

Saved = .true.

NScaleMax = 5
do while ((Qold-Qnew) > CMSThreshold)
  NCMSScale = NCMSScale+1
  if (NCMSScale == nScaleMax) then
    write(u6,'(6X,A)') 'Scaling does not save Qaa from decreasing.'
    write(u6,'(6X,A)') 'Q_a-a decreases for this step.'
    Saved = .false.
    exit
  end if
  R(:) = RCopy(:)
  GDState(:) = GDCopy(:)
  DgState(:) = DgCopy(:)
  X(:) = 0.1_wp*X(:)

  call UpDateRotMat(R,DeltaR,X,lRoots,nSPair)
  call RotGD(GDstate,DeltaR,nGD,lRoots,NAC2)
  call RotGD(Dgstate,DeltaR,nGD,lRoots,NAC2)
  call TransposeMat(Dgorbit,Dgstate,nGD,lRoots2,NAC2)
  call TransposeMat(GDorbit,GDstate,nGD,lRoots2,NAC2)
  call CalcDDg(DDg,GDorbit,Dgorbit,nDDg,nGD,lRoots2,NAC2)
  call CalcQaa(Qnew,DDg,lRoots,nDDg)

end do

return

end subroutine CMSScaleX
