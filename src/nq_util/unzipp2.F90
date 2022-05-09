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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
subroutine UnzipP2(P2Unzip,P2MO,nP2Act)

use nq_Info, only: iOff_Ash, mIrrep, NASH, NASHT, NASHT4
use Index_Functions, only: iTri
use Constants, only: One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nP2Act
real(kind=wp) :: P2Unzip(NASHT4), P2MO(nP2Act)
! Input: nP2Act P2MO
! Output: P2Unzip
integer(kind=iwp) :: I, IAct, iIrrep, IJ, IJKL, IOFF1, IOff2, IOff3, J, JAct, jIrrep, K, kAct, kIrrep, KL, L, LAct, lIrrep, &
                     NASHT2, NASHT3
real(kind=wp) :: Fact

if (NASHT4 == 0) return

NASHT2 = NASHT**2
NASHT3 = NASHT2*NASHT

do IIrrep=0,mIrrep-1
  do I=1,NASH(iIrrep)
    IAct = iOff_Ash(iIrrep)+I
    IOff1 = (IAct-1)*NASHT3
    do jIrrep=0,mIrrep-1
      do J=1,NASH(JIrrep)
        JAct = iOff_Ash(JIrrep)+J
        IOff2 = IOff1+(JAct-1)*NASHT2
        IJ = iTri(IAct,JAct)
        do kIrrep=0,mIrrep-1
          do K=1,NASH(KIrrep)
            KAct = IOff_Ash(KIrrep)+K
            IOff3 = IOff2+(KAct-1)*NASHT
            do lIrrep=0,mIrrep-1
              do L=1,NASH(lIrrep)
                LAct = IOff_Ash(LIrrep)+L
                KL = iTri(KAct,LAct)
                IJKL = iTri(ij,kl)
                Fact = Half
                if ((ij >= kl) .and. (kAct == lAct)) Fact = One
                if ((kl >= ij) .and. (iAct == jAct)) Fact = One
                P2Unzip(IOff3+LAct) = P2MO(ijkl)*Fact
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do

return

end subroutine UnzipP2
