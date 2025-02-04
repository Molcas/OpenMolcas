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

function NDXFSM(NSMOB,NSMSX,MXPOBS,NO1PS,NO2PS,NO3PS,NO4PS,IDXSM,ADSXA,SXDXSX,IS12,IS34,IS1234,IPRNT)
! Number of double excitations with total symmetry IDXSM
!
! IS12 (0,1,-1)   : Permutational symmetry between index 1 and 2
! IS34 (0,1,-1)   : Permutational symmetry between index 3 and 3
! IS1234 (0,1,-1) : permutational symmetry between index 12 and 34

use Definitions, only: u6

! General input
integer ADSXA(MXPOBS,2*MXPOBS), SXDXSX(2*MXPOBS,4*MXPOBS)
! Specific input
integer NO1PS(*), NO2PS(*), NO3PS(*), NO4PS(*)

N12 = 0
N34 = 0
MDX = 0
do I12SM=1,NSMSX
  do I1SM=1,NSMOB
    I2SM = ADSXA(I1SM,I12SM)
    if ((IS12 /= 0) .and. (I1SM < I2SM)) goto 190
    if (IS12 == 0) then
      I12NUM = (I1SM-1)*NSMSX+I2SM
    else
      I12NUM = I1SM*(I1SM+1)/2+I2SM
    end if
    if ((IS12 == 0) .or. (I1SM /= I2SM)) then
      N12 = NO1PS(I1SM)*NO2PS(I2SM)
    else if ((IS12 == 1) .and. (I1SM == I2SM)) then
      N12 = NO1PS(I1SM)*(NO1PS(I1SM)+1)/2
    else if ((IS12 == -1) .and. (I1SM == I2SM)) then
      N12 = NO1PS(I1SM)*(NO1PS(I1SM)-1)/2
    end if
    I34SM = SXDXSX(I12SM,IDXSM)
    do I3SM=1,NSMOB
      I4SM = ADSXA(I3SM,I34SM)
      if ((IS34 /= 0) .and. (I3SM < I4SM)) goto 90
      if (IS34 == 0) then
        I34NUM = (I3SM-1)*NSMSX+I4SM
      else
        I34NUM = I3SM*(I3SM+1)/2+I4SM
      end if
      if ((IS1234 /= 0) .and. (I12NUM < I34NUM)) goto 90
      if ((IS34 == 0) .or. (I3SM /= I4SM)) then
        N34 = NO3PS(I3SM)*NO4PS(I4SM)
      else if ((IS34 == 1) .and. (I3SM == I4SM)) then
        N34 = NO3PS(I3SM)*(NO3PS(I3SM)+1)/2
      else if ((IS34 == -1) .and. (I3SM == I4SM)) then
        N34 = NO3PS(I3SM)*(NO3PS(I3SM)-1)/2
      end if
      if ((IS1234 == 0) .or. (I12NUM /= I34NUM)) then
        MDX = MDX+N12*N34
      else if ((IS1234 == 1) .and. (I12NUM == I34NUM)) then
        MDX = MDX+N12*(N12+1)/2
      else if ((IS1234 == -1) .and. (I12NUM == I34NUM)) then
        MDX = MDX+N12*(N12-1)/2
      end if
      !write(u6,*) ' I1SM I2SM I3SM I4SM MDX'
      !write(u6,*) I1SM,I2SM,I3SM,I4SM,MDX
90  continue
    end do
190 continue
  end do
end do

NDXFSM = MDX

NTEST = 0
NTEST = max(NTEST,IPRNT)
if (NTEST /= 0) write(u6,*) ' Number of double excitations obtained ',MDX

end function NDXFSM
