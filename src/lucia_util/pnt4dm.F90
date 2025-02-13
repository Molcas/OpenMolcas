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

subroutine PNT4DM(NSMOB,MXPOBS,NO1PS,NO2PS,NO3PS,NO4PS,IDXSM,ADSXA,SXDXSX,IS12,IS34,IS1234,IPNTR,ISM4A,ADASX)
! Pointer for 4 dimensionl array with total symmetry IDXSM
! Pointer is given as 3 dimensional array corresponding
! to the first 3 indices
! Symmetry of last index is give by ISM4
!
! IS12 (0,1,-1)   : Permutational symmetry between indices 1 and 2
! IS34 (0,1,-1)   : Permutational symmetry between indices 3 and 3
! IS1234 (0,1,-1) : permutational symmetry between indices 12 and 34

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSMOB, MXPOBS, NO1PS(NSMOB), NO2PS(NSMOB), NO3PS(NSMOB), NO4PS(NSMOB), IDXSM, &
                                 ADSXA(MXPOBS,2*MXPOBS), SXDXSX(2*MXPOBS,4*MXPOBS), IS12, IS34, IS1234, ADASX(MXPOBS,MXPOBS)
integer(kind=iwp), intent(out) :: IPNTR(NSMOB,NSMOB,NSMOB), ISM4A(NSMOB,NSMOB,NSMOB)
integer(kind=iwp) :: I12NUM, I12SM, I1SM, I2SM, I34NUM, I34SM, I3SM, I4SM, IOFF, N12, N34, NTEST

IPNTR(:,:,:) = 0
ISM4A(:,:,:) = 0

!write(u6,*) 'NO1PS NO2PS NO3PS NO4PS'
!call IWRTMA(NO1PS,1,NSMOB,1,NSMOB)
!call IWRTMA(NO2PS,1,NSMOB,1,NSMOB)
!call IWRTMA(NO3PS,1,NSMOB,1,NSMOB)
!call IWRTMA(NO4PS,1,NSMOB,1,NSMOB)
IOFF = 1
N12 = 0
N34 = 0

do I1SM=1,NSMOB
  do I2SM=1,NSMOB
    I12SM = ADASX(I1SM,I2SM)
    I34SM = SXDXSX(I12SM,IDXSM)
    if (I34SM == 0) cycle
    if ((IS12 /= 0) .and. (I1SM < I2SM)) cycle
    if (IS12 == 0) then
      I12NUM = (I1SM-1)*NSMOB+I2SM
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
    do I3SM=1,NSMOB
      I4SM = ADSXA(I3SM,I34SM)
      if (I4SM == 0) cycle
      if ((IS34 /= 0) .and. (I3SM < I4SM)) cycle
      if (IS34 == 0) then
        I34NUM = (I3SM-1)*NSMOB+I4SM
      else
        I34NUM = I3SM*(I3SM+1)/2+I4SM
      end if
      if ((IS1234 /= 0) .and. (I12NUM < I34NUM)) cycle
      if ((IS34 == 0) .or. (I3SM /= I4SM)) then
        N34 = NO3PS(I3SM)*NO4PS(I4SM)
      else if ((IS34 == 1) .and. (I3SM == I4SM)) then
        N34 = NO3PS(I3SM)*(NO3PS(I3SM)+1)/2
      else if ((IS34 == -1) .and. (I3SM == I4SM)) then
        N34 = NO3PS(I3SM)*(NO3PS(I3SM)-1)/2
      end if
      if ((IS1234 == 0) .or. (I12NUM /= I34NUM)) then
        IPNTR(I1SM,I2SM,I3SM) = IOFF
        ISM4A(I1SM,I2SM,I3SM) = I4SM
        IOFF = IOFF+N12*N34
      else if ((IS1234 == 1) .and. (I12NUM == I34NUM)) then
        IPNTR(I1SM,I2SM,I3SM) = IOFF
        ISM4A(I1SM,I2SM,I3SM) = I4SM
        IOFF = IOFF+N12*(N12+1)/2
      else if ((IS1234 == -1) .and. (I12NUM == I34NUM)) then
        IPNTR(I1SM,I2SM,I3SM) = IOFF
        ISM4A(I1SM,I2SM,I3SM) = I4SM
        IOFF = IOFF+N12*(N12-1)/2
      end if
      !write(u6,*) ' I1SM I2SM I3SM I4SM    IOFF'
      !write(u6,'(1X,4I4,I9)') I1SM,I2SM,I3SM,I4SM,IOFF
    end do
  end do
end do

!write(u6,*) ' PNT4DM, 64 elemets of IPNTR'
!call IWRTMA(IPNTR,1,64,1,64)
NTEST = 0
if (NTEST /= 0) write(u6,*) ' Length of 4 index array ',IOFF-1

end subroutine PNT4DM
