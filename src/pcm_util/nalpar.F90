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

function NAlPar(IAt,IAn,NBond,IBond,Chg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NAlPar
integer(kind=iwp), parameter :: MxBond = 12
integer(kind=iwp), intent(in) :: IAt, IAn(*), NBond(*), IBond(MxBond,*)
real(kind=wp), intent(in) :: Chg(*)
integer(kind=iwp) :: ianj, iank, jat, jj, kat, kk, NArI, nbj, nbk, NCSP2J, NHetJ
real(kind=wp) :: chgk

NAlPar = -1
NArI = 0
do jj=1,3
  NCSP2J = 0
  NHetJ = 0
  jat = IBond(jj,iat)
  ianj = ian(jat)
  nbj = nbond(jat)
  if ((ianj == 7) .and. (nbj >= 3)) NCSP2J = 2
  if ((ianj == 6) .and. (nbj == 3)) then
    do kk=1,3
      kat = IBond(kk,JAt)
      iank = ian(KAt)
      nbk = nbond(Kat)
      chgk = chg(KAt)
      if (chgk < 0.4_wp) then
        if ((iank == 6) .and. (nbk == 3)) NCSP2J = NCSP2J+1
        if ((iank == 8) .or. (iank == 9)) NHetJ = NHetJ+1
        if ((iank == 17) .or. (iank == 35) .or. (iank == 53)) NHetJ = NHetJ+1
        if (iank == 7) then
          if (nbk <= 2) NHetJ = NHetJ+1
          if (nbk >= 3) NCSP2J = NCSP2J+1
        end if
      end if
    end do
  end if
  if ((NCSP2J >= 2) .and. (NHetJ == 0)) NarI = NArI+1
end do
if (NArI >= 2) NAlPar = 1

return

end function NAlPar
