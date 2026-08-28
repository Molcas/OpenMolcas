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

subroutine GETINC_ABT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,INTLST,ICOUL,ieaw)
! Obtain integrals
! ICOUL = 0 :      XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
!                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
! ICOUL = 1 :      XINT(IJ,KL) = (IJ!KL)
!
! Version for integrals stored in INTLST

use Index_Functions, only: iTri
use MCLR_Data, only: IBTSOB, NTSOB
use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: XINT(*)
integer(kind=iwp), intent(in) :: ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IXCHNG, IKSM, JLSM, ICOUL, ieaw
real(kind=wp), intent(in) :: Intlst(*)
integer(kind=iwp) :: iBas, iInt, IJ, IJKL, IL, IMIN, iOff, iOrb, jBas, JK, JMIN, jOff, jOrb, kBas, KL, kOff, kOrb, lBas, lOff, lOrb
real(kind=wp) :: SGN

iOrb = NTSOB(ITP,ISM)
jOrb = NTSOB(JTP,JSM)
kOrb = NTSOB(KTP,KSM)
lOrb = NTSOB(LTP,LSM)
iOff = IBTSOB(ITP,ISM)
jOff = IBTSOB(JTP,JSM)
kOff = IBTSOB(KTP,KSM)
lOff = IBTSOB(LTP,LSM)

if (ICOUL == 0) then

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
          ij = iTri(iBas,jBas)
          kl = iTri(kBas,lBas)
          SGN = One
          if ((ij < kl) .and. (ieaw /= 0)) SGN = -One
          ijkl = iTri(ij,kl)
          Xint(iInt) = SGN*Intlst(ijkl)
          iInt = iInt+1
        end do
      end do
    end do
  end do

  ! Collect Exchange terms

  if (IXCHNG /= 0) then
    iint = 1
    do lBas=lOff,lOff+lOrb-1
      jMin = jOff
      if (JLSM /= 0) jMin = lBas
      do jBas=jMin,jOff+jOrb-1
        do kBas=kOff,kOff+kOrb-1
          iMin = iOff
          if (IKSM /= 0) iMin = kBas
          do iBas=iMin,iOff+iOrb-1
            il = iTri(iBas,lBas)
            jk = iTri(jBas,kBas)
            ijkl = iTri(il,jk)
            XInt(iInt) = XInt(iInt)-Intlst(ijkl)
            iInt = iInt+1
          end do
        end do
      end do
    end do
  end if
else
  iint = 1
  do lBas=lOff,lOff+lOrb-1
    do kBas=kOff,kOff+kOrb-1
      do jBas=joff,jOff+jOrb-1
        do iBas=ioff,iOff+iOrb-1
          ij = iTri(iBas,jBas)
          kl = iTri(kBas,lBas)
          ijkl = iTri(ij,kl)
          Xint(iInt) = Intlst(ijkl)
          iInt = iint+1
        end do
      end do
    end do
  end do
end if

end subroutine GETINC_ABT
