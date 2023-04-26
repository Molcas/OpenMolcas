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

subroutine DefParo2v4(NaGrp,NbeGrp,NaSGrp,NbeSgrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
! This routine does:
! define parameters using NaGrp,NbeGrp,NaSGrp,NbeSgrp
!
! I/O parameter description:
! NaGrp    - # of groups in a set (I)
! NbeGrp   - # of groups in be set (I)
! NaSGrp   - # of subgroups in each (a)' group (I)
! NbeSGrp  - # of subgroups in each (be)' group (I)
! mdGrpa   - # maximal dimension od (a)' group (O)
! mdGrpbe  - # maximal dimension od (be)' group (O)
! mdSGrpa  - # maximal dimension od (a)" subgroup (O)
! mdSGrpbe - # maximal dimension od (be)" subgroup (O)

use chcc_global, only: DimGrpa, DimGrpbe, DimSGrpa, DimSGrpbe, GrpaLow, GrpaUp, GrpbeLow, GrpbeUp, L2Name, maxGrp, maxSGrp, nv, &
                       Tmp3Name
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NaGrp, NbeGrp, NaSGrp, NbeSgrp
integer(kind=iwp), intent(out) :: mdGrpa, mdGrpbe, mdSGrpa, mdSGrpbe
integer(kind=iwp) :: i, ij, j, LowGrpa(maxGrp), LowGrpbe(maxGrp), LowSGrpa(maxSGrp), LowSGrpbe(maxSGrp), UpGrpa(maxGrp), &
                     UpGrpbe(maxGrp), UpSGrpa(maxSGrp), UpSGrpbe(maxSGrp)
real(kind=wp) :: rdim

!1 define parameters of Groups of a(b) set

rdim = real(nv,kind=wp)/real(NaGrp,kind=wp)

Grpalow(1) = 1
Grpaup(1) = NaSGrp

do i=1,NaGrp

  if (i == 1) then
    UpGrpa(i) = int(rdim*i)
    LowGrpa(i) = 1
  else if (i == NaGrp) then
    UpGrpa(i) = nv
    LowGrpa(i) = UpGrpa(i-1)+1
  else
    UpGrpa(i) = int(rdim*i)
    LowGrpa(i) = UpGrpa(i-1)+1
  end if

  DimGrpa(i) = (UpGrpa(i)-LowGrpa(i))+1

  if (i > 1) then
    Grpalow(i) = Grpalow(i-1)+NaSGrp
    Grpaup(i) = Grpaup(i-1)+NaSGrp
  end if

end do

!2 define parameters of Groups of be(ga) set

rdim = real(nv,kind=wp)/real(NbeGrp,kind=wp)

Grpbelow(1) = 1
Grpbeup(1) = NbeSGrp

do i=1,NbeGrp

  if (i == 1) then
    UpGrpbe(i) = int(rdim*i)
    LowGrpbe(i) = 1
  else if (i == NbeGrp) then
    UpGrpbe(i) = nv
    LowGrpbe(i) = UpGrpbe(i-1)+1
  else
    UpGrpbe(i) = int(rdim*i)
    LowGrpbe(i) = UpGrpbe(i-1)+1
  end if

  DimGrpbe(i) = (UpGrpbe(i)-LowGrpbe(i))+1

  if (i > 1) then
    Grpbelow(i) = Grpbelow(i-1)+NbeSGrp
    Grpbeup(i) = Grpbeup(i-1)+NbeSGrp
  end if

end do

!3 define parameters of SubGroups of a(b) set

ij = 0
do j=1,NaGrp

  rdim = real(DimGrpa(j),kind=wp)/real(NaSGrp,kind=wp)

  do i=1,NaSGrp
    ij = ij+1

    if (i == 1) then
      UpSGrpa(ij) = (LowGrpa(j)-1)+int(rdim*i)
      LowSGrpa(ij) = LowGrpa(j)
    else if (i == NaSGrp) then
      UpSGrpa(ij) = UpGrpa(j)
      LowSGrpa(ij) = UpSGrpa(ij-1)+1
    else
      UpSGrpa(ij) = (LowGrpa(j)-1)+int(rdim*i)
      LowSGrpa(ij) = UpSGrpa(ij-1)+1
    end if

    DimSGrpa(ij) = (UpSGrpa(ij)-LowSGrpa(ij))+1

  end do
end do

!4 define parameters of SubGroups of be(ga) set

ij = 0
do j=1,NbeGrp

  rdim = real(DimGrpbe(j),kind=wp)/real(NbeSGrp,kind=wp)

  do i=1,NbeSGrp
    ij = ij+1

    if (i == 1) then
      UpSGrpbe(ij) = (LowGrpbe(j)-1)+int(rdim*i)
      LowSGrpbe(ij) = LowGrpbe(j)
    else if (i == NbeSGrp) then
      UpSGrpbe(ij) = UpGrpbe(j)
      LowSGrpbe(ij) = UpSGrpbe(ij-1)+1
    else
      UpSGrpbe(ij) = (LowGrpbe(j)-1)+int(rdim*i)
      LowSGrpbe(ij) = UpSGrpbe(ij-1)+1
    end if

    DimSGrpbe(ij) = (UpSGrpbe(ij)-LowSGrpbe(ij))+1

  end do
end do

!5 find maximal dimensions mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe

mdGrpa = DimGrpa(1)
do i=1,NaGrp
  if (DimGrpa(i) > mdGrpa) mdGrpa = DimGrpa(i)
end do

mdGrpbe = DimGrpbe(1)
do i=1,NbeGrp
  if (DimGrpbe(i) > mdGrpbe) mdGrpbe = DimGrpbe(i)
end do

mdSGrpa = DimSGrpa(1)
do ij=1,NaGrp*NaSGrp
  if (DimSGrpa(ij) > mdSGrpa) mdSGrpa = DimSGrpa(ij)
end do

mdSGrpbe = DimSGrpbe(1)
do ij=1,NbeGrp*NbeSGrp
  if (DimSGrpbe(ij) > mdSGrpbe) mdSGrpbe = DimSGrpbe(ij)
end do

!6 def L2Names

do i=1,NaGrp
  do j=1,NbeGrp
    call DefParo3v3Hlp1(i,j,'L2',L2Name(i,j))
  end do
end do

!7 def Tmp3Names

do i=1,MaxSGrp
  do j=1,MaxSGrp
    call DefParo3v3Hlp1(i,j,'X3',Tmp3Name(i,j))
  end do
end do

return

end subroutine DefParo2v4
