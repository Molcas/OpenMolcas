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

subroutine distnodes()
! this routine distributes nodes to different parts

use ccsd_global, only: idaaaa, idaabb, idab, idabba, idbaab, idbbaa, idbbbb, ideffab, idfin, nprocab
use Para_Info, only: nProcs
use Constants, only: One, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i
real(kind=wp) :: efftot

!tmp ta zatial takto
if (nProcs == 1) then

  !I def nodes for sumoverab
  nprocab = 1
  idab(1) = 0
  ideffab(1) = One

  !II def nodes for sumoverb and intermezzo
  idaaaa = 0
  idbaab = 0
  idbbaa = 0
  idbbbb = 0
  idaabb = 0
  idabba = 0

  !III def node for finale
  idfin = 0

else if (nProcs == 2) then

  !I def nodes for sumoverab
  nprocab = 1
  idab(1) = 0
  ideffab(1) = One

  !II def nodes for sumoverb and intermezzo
  idaaaa = 1
  idbaab = 1
  idbbaa = 1
  idbbbb = 1
  idaabb = 1
  idabba = 1

  !III def node for finale
  idfin = 1

else if (nProcs == 3) then

  !I def nodes for sumoverab
  nprocab = 1
  idab(1) = 0
  ideffab(1) = One

  !II def nodes for sumoverb and intermezzo
  idaaaa = 1
  idbaab = 1
  idbbaa = 1
  idbbbb = 2
  idaabb = 2
  idabba = 2

  !III def node for finale
  idfin = 1

else if (nProcs == 4) then

  !I def nodes for sumoverab
  nprocab = 4
  idab(1) = 0
  idab(2) = 1
  idab(3) = 2
  idab(4) = 3
  ideffab(1) = Quart
  ideffab(2) = Quart
  ideffab(3) = Quart
  ideffab(4) = Quart

  !II def nodes for sumoverb and intermezzo
  idaaaa = 0
  idbaab = 1
  idbbaa = 1
  idbbbb = 2
  idaabb = 2
  idabba = 3

  !III def node for finale
  idfin = 3

else if (nProcs == 5) then

  !I def nodes for sumoverab
  nprocab = 1
  idab(1) = 0
  ideffab(1) = One

  !II def nodes for sumoverb and intermezzo
  idaaaa = 1
  idbaab = 1
  idbbaa = 2
  idbbbb = 3
  idaabb = 3
  idabba = 4

  !III def node for finale
  idfin = 2

else if (nProcs == 6) then

  !I def nodes for sumoverab
  nprocab = 6
  idab(1) = 0
  idab(2) = 1
  idab(3) = 2
  idab(4) = 3
  idab(5) = 4
  idab(6) = 5
  ideffab(1) = One
  ideffab(2) = One
  ideffab(3) = One
  ideffab(4) = One
  ideffab(5) = One
  ideffab(6) = One

  !II def nodes for sumoverb and intermezzo
  idaaaa = 0
  idbaab = 1
  idbbaa = 2
  idbbbb = 3
  idaabb = 4
  idabba = 5

  !III def node for finale
  idfin = 3

else if (nProcs == 10) then

  !I def nodes for sumoverab
  nprocab = 4
  idab(1) = 0
  idab(2) = 1
  idab(3) = 2
  idab(4) = 3
  ideffab(1) = One
  ideffab(2) = One
  ideffab(3) = One
  ideffab(4) = One

  !II def nodes for sumoverb and intermezzo
  idaaaa = 4
  idbaab = 5
  idbbaa = 6
  idbbbb = 7
  idaabb = 8
  idabba = 9

  !III def node for finale
  idfin = 5

else

  !I def nodes for sumoverab
  nprocab = nProcs
  do i=1,nprocab
    idab(i) = i-1
  end do
  ideffab(1:nprocab) = One

  !II def nodes for sumoverb and intermezzo
  idaaaa = 0
  idbaab = 1
  idbbaa = 2
  idbbbb = 3
  idaabb = 4
  idabba = 5

  !III def node for finale
  idfin = 6

end if

return

!tmp koniec tmp riesenia

if (nProcs == 1) then

  !I def nodes for sumoverab
  nprocab = 1
  idab(1) = 0
  ideffab(1) = One

  !II def nodes for sumoverb and intermezzo
  idaaaa = 0
  idbaab = 0
  idbbaa = 0
  idbbbb = 0
  idaabb = 0
  idabba = 0

  !III def node for finale
  idfin = 0

else if (nProcs == 2) then

  !I def nodes for sumoverab
  nprocab = 2
  idab(1) = 0
  idab(2) = 1
  ideffab(1) = Half
  ideffab(2) = Half

  !II def nodes for sumoverb and intermezzo
  idaaaa = 0
  idbaab = 0
  idbbaa = 0
  idbbbb = 1
  idaabb = 1
  idabba = 1

  !III def node for finale
  idfin = 0

else if (nProcs == 3) then

  !I def nodes for sumoverab
  nprocab = 3
  idab(1) = 0
  idab(2) = 1
  idab(3) = 2
  ideffab(1) = 0.333_wp
  ideffab(2) = 0.333_wp
  ideffab(3) = 0.333_wp

  !II def nodes for sumoverb and intermezzo
  idaaaa = 1
  idbaab = 1
  idbbaa = 1
  idbbbb = 2
  idaabb = 2
  idabba = 2

  !III def node for finale
  idfin = 0

else if (nProcs == 4) then

  !I def nodes for sumoverab
  nprocab = 2
  idab(1) = 2
  idab(2) = 3
  ideffab(1) = One
  ideffab(2) = One

  !II def nodes for sumoverb and intermezzo
  idaaaa = 0
  idbaab = 0
  idbbaa = 0
  idbbbb = 1
  idaabb = 1
  idabba = 1

  !III def node for finale
  idfin = 0

else

  !I def nodes for sumoverab
  nprocab = nProcs-2
  idab(1) = 1
  idab(2) = 2
  idab(3) = 4
  idab(4) = 5
  idab(5) = 6
  idab(6) = 7
  ideffab(1:nprocab) = One/nprocab
  ideffab(1) = Half*ideffab(1)
  ideffab(2) = Half*ideffab(2)

  !II def nodes for sumoverb and intermezzo
  idaaaa = 1
  idbaab = 0
  idbbaa = 0
  idbbbb = 2
  idaabb = 3
  idabba = 3

  !III def node for finale
  idfin = 1

end if

! renormalize ideffab

efftot = sum(ideffab(1:nprocab))

ideffab(1:nprocab) = ideffab(1:nprocab)/efftot

return

end subroutine distnodes
