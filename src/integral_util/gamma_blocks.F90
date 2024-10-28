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

subroutine Gamma_Blocks(iTable,nBlocks,nIrrep)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nBlocks, nIrrep
integer(kind=iwp), intent(out) :: iTable(6,nBlocks)
integer(kind=iwp) :: iBlock, IBOT, idid(8), IND, IRREP, IRREP1, IRREP2, IRREP3, IRREP4, ITMP

iBlock = 0

! Observe the loop ordering straight from GAMSRT in Aces 2.
! THIS LOOP ORDER DOES NOT COMPLY WITH THE MANUAL!

! AAAA

!write(u6,*) 'AAAA'
IND = 0
do IRREP=1,NIRREP
  IND = IND+1
  iBlock = iBlock+1
  iTable(1,iBlock) = 1
  iTable(2:5,iBlock) = IRREP-1
  iTable(6,iBlock) = IND
  !write(u6,*) IND,IRREP-1,IRREP-1,IRREP-1,IRREP-1
end do

! AABB

!write(u6,*) 'AABB'
IND = 0
do IRREP1=2,NIRREP
  do IRREP2=1,IRREP1-1
    IND = IND+1
    iBlock = iBlock+1
    iTable(1,iBlock) = 2
    iTable(2:3,iBlock) = IRREP2-1
    iTable(4:5,iBlock) = IRREP1-1
    iTable(6,iBlock) = IND
    !write(u6,*) IND,IRREP2-1,IRREP2-1,IRREP1-1,IRREP1-1
  end do
end do

! ABAB

!write(u6,*) 'ABAB'
IND = 0
do IRREP=2,NIRREP
  do IRREP1=1,NIRREP
    IRREP2 = ieor(IRREP-1,IRREP1-1)+1
    if (IRREP2 > IRREP1) then
      IND = IND+1
      iBlock = iBlock+1
      iTable(1,iBlock) = 3
      iTable(2:4:2,iBlock) = IRREP1-1
      iTable(3:5:2,iBlock) = IRREP2-1
      iTable(6,iBlock) = IND
      !write(u6,*) IND,IRREP1-1,IRREP2-1,IRREP1-1,IRREP2-1
    end if
  end do
end do

! ABCD

!write(u6,*) 'ABCD'
IND = 0
do IRREP=2,NIRREP
  do IRREP1=1,NIRREP
    IRREP2 = ieor(IRREP-1,IRREP1-1)+1
    if (IRREP2 < IRREP1) cycle
    IBOT = max(IRREP1,IRREP2)+1
    IDID(:) = 0
    do ITMP=IBOT,NIRREP
      IRREP4 = ieor(ITMP-1,IRREP-1)+1
      IRREP3 = min(ITMP,IRREP4)
      IRREP4 = max(ITMP,IRREP4)
      if (max(IDID(IRREP4),IDID(IRREP3)) /= 0) cycle
      IDID(IRREP3) = 1
      IDID(IRREP4) = 1
      IND = IND+1
      iBlock = iBlock+1
      iTable(1,iBlock) = 4
      iTable(2,iBlock) = IRREP3-1
      iTable(3,iBlock) = IRREP4-1
      iTable(4,iBlock) = IRREP1-1
      iTable(5,iBlock) = IRREP2-1
      iTable(6,iBlock) = IND
      !write(u6,*) IND,IRREP3-1,IRREP4-1,IRREP1-1,IRREP2-1
    end do
  end do
end do

end subroutine Gamma_Blocks
