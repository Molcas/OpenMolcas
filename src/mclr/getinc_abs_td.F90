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

subroutine GETINC_ABS_td(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IKSM,JLSM,INTLST,ICTL)
! Obtain integrals
! ICOUL = 0 :      XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
!                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
! ICOUL = 1 :      XINT(IJ,KL) = (IJ!KL)
!
! Version for integrals stored in INTLST

use Index_Functions, only: iTri
use MCLR_Data, only: IBTSOB, NACOB, NTSOB
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: XINT(*)
integer(kind=iwp), intent(in) :: ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IKSM, JLSM, ICTL
real(kind=wp), intent(in) :: Intlst(*)
integer(kind=iwp) :: iBas, iInt, IJ, IJKL, IL, ILJK, IMIN, iOff, iOrb, jBas, JK, JMIN, jOff, jOrb, kBas, KL, kOff, kOrb, lBas, &
                     lOff, lOrb, NTASH

iOrb = NTSOB(ITP,ISM)
jOrb = NTSOB(JTP,JSM)
kOrb = NTSOB(KTP,KSM)
lOrb = NTSOB(LTP,LSM)
iOff = IBTSOB(ITP,ISM)
jOff = IBTSOB(JTP,JSM)
kOff = IBTSOB(KTP,KSM)
lOff = IBTSOB(LTP,LSM)
ntash = nacob

if (ICTL == 1) then

  ! Collect Coulomb terms

  iint = 1
  do lBas=lOff,lOff+lOrb-1
    jMin = jOff
    if (JLSM /= 0) jMin = lBas
    do jBas=jMin,jOff+jOrb-1
      do kBas=kOff,kOff+kOrb-1
        iMin = iOff
        if (IKSM /= 0) iMin = kBas
        do iBas=iMin,iOff+iOrb-1
          ij = jbas+ntash*(ibas-1)
          kl = lbas+ntash*(kbas-1)
          ijkl = iTri(ij,kl)
          Xint(iInt) = Intlst(ijkl)
          iInt = iInt+1
        end do
      end do
    end do
  end do

else if (ICTL == 4) then

  ! Collect Coulomb-Exchange terms

  iint = 1
  do lBas=lOff,lOff+lOrb-1
    jMin = jOff
    if (JLSM /= 0) jMin = lBas
    do jBas=jMin,jOff+jOrb-1
      do kBas=kOff,kOff+kOrb-1
        iMin = iOff
        if (IKSM /= 0) iMin = kBas
        do iBas=iMin,iOff+iOrb-1
          il = ibas+ntash*(lbas-1)
          jk = kbas+ntash*(jbas-1)
          iljk = iTri(il,jk)
          ij = ibas+ntash*(jbas-1)
          kl = kbas+ntash*(lbas-1)
          ijkl = iTri(ij,kl)
          XInt(iInt) = Intlst(ijkl)-Intlst(iljk)
          iInt = iInt+1
        end do
      end do
    end do
  end do

else
  call Abend()
end if

end subroutine GETINC_ABS_td
