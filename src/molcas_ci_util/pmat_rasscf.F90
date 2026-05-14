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

subroutine PMAT_RASSCF(P,X)
! RASSCF version IBM-3090: SX section
!
! Purpose: To construct from a canonically ordered list of
!          2-matrix elements a list ordered as the transformed
!          two-electron integrals. P is the input matrix and X is
!          the output matrix. the matrix is multiplied by two and
!          used to construct the Q-matrix in fock.
!
! ********** IBM-3090 MOLCAS Release: 90 02 22 **********

use Symmetry_Info, only: Mul
use Index_Functions, only: iTri, nTri_Elem
use rasscf_global, only: ISTORP
use general_data, only: NASH, NSYM
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: P(*)
real(kind=wp), intent(_OUT_) :: X(*)
integer(kind=iwp) :: IAT, IAU, IAV, IAX, INDF, INDX, LAT, LAU, LAV, LAX, LPMAT, LROW, NAP, NAQ, NAR, NAS, NAT, NAU, NAV, NAX, &
                     NAXE, NSP, NSPQ, NSQ, NSR, NSS, NSS1, NSSM, NTU, NTUVX, NUVX, NVX
real(kind=wp) :: FAC

#ifdef _DEBUGPRINT_
write(u6,*) ' Entering PMAT_RASSCF'
#endif

! Loop over all reordered 2-matrix elements.

LPMAT = ISTORP(NSYM+1)
X(1:LPMAT) = Zero

IAT = 0
do NSP=1,NSYM
  NAP = NASH(NSP)
  if (NAP == 0) cycle
  INDF = ISTORP(NSP)
  NUVX = (ISTORP(NSP+1)-INDF)/NAP
  LROW = 0
  IAU = 0
  do NSQ=1,NSYM
    NAQ = NASH(NSQ)
    if (NAQ == 0) cycle
    NSPQ = Mul(NSP,NSQ)
    IAV = 0
    do NSR=1,NSYM
      NAR = NASH(NSR)
      if (NAR == 0) cycle
      NSS = Mul(NSPQ,NSR)
      if (NSS <= NSR) then
        NAS = NASH(NSS)
        if (NAS /= 0) then
          IAX = 0
          if (NSS /= 1) then
            NSSM = NSS-1
            do NSS1=1,NSSM
              IAX = IAX+NASH(NSS1)
            end do
          end if
          do NAV=1,NAR
            LAV = NAV+IAV
            NAXE = NAS
            if (NSR == NSS) NAXE = NAV
            do NAX=1,NAXE
              LAX = NAX+IAX
              do NAU=1,NAQ
                LAU = NAU+IAU
                LROW = LROW+1
                INDX = INDF+LROW-NUVX
                do NAT=1,NAP
                  INDX = INDX+NUVX
                  LAT = NAT+IAT

                  ! Compute canonical index ntuvx and find prefactor

                  NTU = iTri(LAT,LAU)
                  NVX = nTri_Elem(LAV-1)+LAX
                  NTUVX = iTri(NTU,NVX)
                  FAC = Two
                  if (NTU < NVX) then
                    if ((LAT == LAU) .and. (LAV /= LAX)) then
                      FAC = Four
                    else if ((LAT /= LAU) .and. (LAV == LAX)) then
                      FAC = One
                    end if
                  end if
                  X(INDX) = FAC*P(NTUVX)
                end do
              end do
            end do
          end do
        end if
      end if
      IAV = IAV+NAR
    end do
    IAU = IAU+NAQ
  end do
  IAT = IAT+NAP
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Reordered 2-matrix:'
write(u6,'(1X,10F10.6)') X(1:LPMAT)
#endif

end subroutine PMAT_RASSCF
