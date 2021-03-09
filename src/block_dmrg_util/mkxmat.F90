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
! Copyright (C) 2014, Naoki Nakatani                                   *
!***********************************************************************

subroutine MKXMAT(TORB,XMAT)
! Make full transformation matrix for active space
! from that stored for each symmetry in IAD1M(4)
! Written by N. Nakatani, Oct. 2014

use Definitions, only: wp, iwp

implicit none
#include "rasdim.fh"
#include "caspt2.fh"
real(kind=wp), intent(in) :: TORB(NTORB)
real(kind=wp), intent(out) :: XMAT(NASHT,NASHT)
integer(kind=iwp) :: I, IR1, IR2, IR3, ISTART, ISYM, ITO, ITOEND, ITOSTA, J, JR1, JR2, JR3, NA, NI, NR1, NR2, NR3, NS

if (NASHT > 0) then
  ITOEND = 0
  do ISYM=1,NSYM
    NI = NISH(ISYM)
    NA = NASH(ISYM)
    NR1 = NRAS1(ISYM)
    NR2 = NRAS2(ISYM)
    NR3 = NRAS3(ISYM)
    NS = NSSH(ISYM)
    ITOSTA = ITOEND+1
    ITOEND = ITOEND+NI**2+NR1**2+NR2**2+NR3**2+NS**2

    ! Normally, NRAS2 should only be non-zero, but NRAS1 and NRAS3 should be zero
    ! for DMRG-CASSCF calculation
    ITO = ITOSTA+NI**2
    !ITO = ITOSTA+NI**2-1
    if (NA > 0) then
      ! RAS1
      if (NR1 > 0) then
        ISTART = NAES(ISYM)
        do JR1=1,NR1
          J = ISTART+JR1
          do IR1=1,NR1
            I = ISTART+IR1
            XMAT(I,J) = TORB(ITO)
            ITO = ITO+1
          end do
        end do
      end if
      ! RAS2
      if (NR2 > 0) then
        ISTART = NAES(ISYM)+NR1
        do JR2=1,NR2
          J = ISTART+JR2
          do IR2=1,NR2
            I = ISTART+IR2
            XMAT(I,J) = TORB(ITO)
            ITO = ITO+1
          end do
        end do
      end if
      ! RAS3
      if (NR3 > 0) then
        ISTART = NAES(ISYM)+NR1+NR2
        do JR3=1,NR3
          J = ISTART+JR3
          do IR3=1,NR3
            I = ISTART+IR3
            XMAT(I,J) = TORB(ITO)
            ITO = ITO+1
          end do
        end do
      end if
    end if
  end do
end if

return

end subroutine MKXMAT
