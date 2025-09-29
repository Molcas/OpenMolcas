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

subroutine GETINC_ABS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,INTLST,ICOUL)
! Obtain integrals
! ICOUL = 0 :      XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
!                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
! ICOUL = 1 :      XINT(IJ,KL) = (IJ!KL)
!
! Version for integrals stored in INTLST

use MCLR_Data, only: IBTSOB, NACOB, NTSOB
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: XINT(*)
integer(kind=iwp), intent(in) :: ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IXCHNG, ICOUL
real(kind=wp), intent(in) :: Intlst(*)
integer(kind=iwp) :: iInt, iOff, iOrb, jBas, jInt, jOff, jOrb, kBas, kOff, kOrb, lBas, lOff, lOrb

iOrb = NTSOB(ITP,ISM)
jOrb = NTSOB(JTP,JSM)
kOrb = NTSOB(KTP,KSM)
lOrb = NTSOB(LTP,LSM)
iOff = IBTSOB(ITP,ISM)
jOff = IBTSOB(JTP,JSM)
kOff = IBTSOB(KTP,KSM)
lOff = IBTSOB(LTP,LSM)

! Collect Coulomb terms

if (ICOUL == 0) then
  iInt = 1
  do lBas=lOff,lOff+lOrb-1
    do jBas=jOff,jOff+jOrb-1
      do kBas=kOff,kOff+kOrb-1
        jInt = (lBas-1)*nACOB**3+(kBas-1)*nACOB**2+(jBas-1)*nACOB
        Xint(iInt+iOff:iInt+iOff+iOrb-1) = Intlst(jInt+iOff:jInt+iOff+iOrb-1)
        iInt = iInt+iOrb
      end do
    end do
  end do

  ! Collect Exchange terms

  if (IXCHNG /= 0) then
    iInt = 1
    do lBas=lOff,lOff+lOrb-1
      do jBas=jOff,jOff+jOrb-1
        do kBas=kOff,kOff+kOrb-1
          jInt = (jBas-1)*nACOB**3+(kBas-1)*nACOB**2+(lBas-1)*nACOB
          XInt(iInt+iOff:iInt+iOff+iOrb-1) = XInt(iInt+iOff:iInt+iOff+iOrb-1)-Intlst(jInt+iOff:jInt+iOff+iOrb-1)
          iInt = iInt+iOrb
        end do
      end do
    end do
  end if
else
  iInt = 0
  do lBas=lOff,lOff+lOrb-1
    do kBas=kOff,kOff+kOrb-1
      do jBas=jOff,jOff+jOrb-1
        jInt = (lBas-1)*nACOB**3+(kBas-1)*nACOB**2+(jBas-1)*nACOB
        Xint(iInt+iOff:iInt+iOff+iOrb-1) = Intlst(jInt+iOff:jInt+iOff+iOrb-1)
        iInt = iInt+iOrb
      end do
    end do
  end do
end if

end subroutine GETINC_ABS
