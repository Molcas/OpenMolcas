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

function NCAlph(IAt,NHI,NCSP3I,IAn,NBond,IBond,Chg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NCAlph
integer(kind=iwp), parameter :: MxBond = 12
integer(kind=iwp), intent(in) :: IAt, NHI, NCSP3I, IAn(*), NBond(*), IBond(MxBond,*)
real(kind=wp), intent(in) :: Chg(*)
integer(kind=iwp) :: ianj, iank, iplj, jat, jj, kat, kk, nbj, nbk, NCM1, NCP1, ncsp3j, NHetI, nhetj, nhj

NHetI = 4-NHI-NCSP3I
NCP1 = 0
NCM1 = 0
!write(IOut,'("Atom",I2," NHet=",I1)') IAt,NHetI
do jj=1,4
  jat = IBond(jj,iat)
  ianj = ian(jat)
  nbj = nbond(jat)
  nhj = 0
  ncsp3j = 0
  nhetj = 0
  iplj = 0
  if ((ianj == 6) .and. (nbj == 4)) then
    do kk=1,4
      kat = IBond(kk,JAt)
      iank = ian(KAt)
      nbk = nbond(kat)
      if (iank == 1) NHJ = NHJ+1
      if ((iank == 6) .and. (nbk == 4)) NCSP3J = NCSP3J+1
      if (chg(kat) > 0.4_wp) iplj = 1
    end do
    NHETJ = 4-NHJ-NCSP3J
    if ((NHetI >= 0) .and. (NHetJ == 0)) NCP1 = NCP1+1
    if ((NHetI == 0) .and. (NHetJ > 0) .and. (iplj == 0)) NCM1 = NCM1+1
  end if
  !write(IOut,*) JAt,NHetJ,NCP1,NCM1
end do
NCAlph = NCP1-NCM1

return

end function NCAlph
