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

subroutine mkTraCI(nTORB,TORB,STSYM,nConf,CI)

use caspt2_module, only: nAES, nAsh, nIsh, nRas1, nRas2, nRas3, nSsh, nSym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTORB, STSYM, nCONF
real(kind=wp), intent(inout) :: TORB(nTORB), CI(nConf)
integer(kind=iwp) :: iSTART, ISYM, ITO, ITOEND, ITOSTA, NA, NI, NR1, NR2, NR3, NS

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

  ITO = ITOSTA+NI**2
  if (NA <= 0) cycle
  if (NR1 > 0) then
    ISTART = NAES(ISYM)+1
    call TRACI_RPT2(ISTART,NR1,TORB(ITO),STSYM,NCONF,CI)
  end if
  ITO = ITO+NR1**2
  if (NR2 > 0) then
    ISTART = NAES(ISYM)+NR1+1
    call TRACI_RPT2(ISTART,NR2,TORB(ITO),STSYM,NCONF,CI)
  end if
  ITO = ITO+NR2**2
  if (NR3 > 0) then
    ISTART = NAES(ISYM)+NR1+NR2+1
    call TRACI_RPT2(ISTART,NR3,TORB(ITO),STSYM,NCONF,CI)
  end if
end do

end subroutine mkTraCI
