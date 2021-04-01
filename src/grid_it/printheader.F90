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

subroutine PrintHeader(nMOs,nShowMOs,nShowMOs_ab,nCoor,nInc,iiCoord,nTypes,iCRSIZE,NBYTES,NINLINE,nBlocks)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use grid_it_globals, only: Cutoff, GridAxis1, GridAxis2, GridAxis3, GridOrigin, iGridNpt, isCutOff, isLine, isLuscus, isMOPack, &
                           isTheOne, isUHF, LID, LID_ab, LuVal, LuVal_ab, VERSION
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nMOs, nShowMOs, nShowMOs_ab, nCoor, nInc, iiCoord, nTypes(7), iCRSIZE, NBYTES, NINLINE
integer(kind=iwp), intent(inout) :: nBlocks
integer(kind=iwp) :: i, iiUHF, LuVal_, NFIRST, NLAST, nn, nn1, nShowMOs_
character(len=128) :: line

LuVal_ = LuVal
if (isLuscus) LuVal_ = LID
nShowMOs_ = nShowMOs
!call bXML('Header')
!call iXML('UHF',merge(1,0,isUHF))
do iiUHF=0,merge(1,0,isUHF)
  if (iiUHF == 1) then
    LuVal_ = LuVal_ab
    if (isLuscus) LuVal_ = LID_ab
    nShowMOs_ = nShowMOs_ab
  end if
  if (isLuscus) then
    ! debug
    !write(u6,*) 'N_of_MO=',nMOs
    !write(u6,*) 'N_of_Grids=',nShowMOs_
    !write(u6,*) 'N_of_Points=',nCoor
    !write(u6,*) 'Block_Size=',nInc
    !write(u6,*) 'N_Blocks=',nBlocks
    !write(u6,*) 'Is_cutoff=',isCutOff
    !write(u6,*) 'CutOff=',Cutoff
    !write(u6,*) 'N_P=',iiCoord
    ! end debug

    if (nMOs > 99999) then
      write(u6,*) 'Number of MO''s can''t be larger that 99999'
      call Quit_OnUserError()
    end if
    if (nShowMOs_ > 9999) then
      write(u6,*) 'Number of grids can''t be larger that 9999'
      call Quit_OnUserError()
    end if
    if (nCoor > 99999999) then
      write(u6,*) 'Number of points can''t be larger that 99999999'
      call Quit_OnUserError()
    end if
    if (nInc > 999999) then
      write(u6,*) 'Block size can''t be larger that 999999'
      call Quit_OnUserError()
    end if
    nBlocks = nCoor/nInc+1
    if (nBlocks > 9999) then
      write(u6,*) 'Number of blocks can''t be larger that 9999'
      !write(u6,*) 'nBlocks = ',nBlocks
      call Quit_OnUserError()
    end if
    !if (isCutOff > 9) then
    !  write(u6,*) 'Wrong cutoff'
    !  call Quit_OnUserError()
    !end if
    if (iiCoord > 99999999) then
      write(u6,*) 'N_P can''t be larger that 99999999'
      call Quit_OnUserError()
    end if
    write(LINE,'(''<GRID>'')')
    call PRINTLINE(LUVAL_,LINE,6,.false.)

    write(LINE,1000) nMOs,nShowMOs_,nCoor,nInc,nBlocks,merge(1,0,isCutOff),Cutoff,iiCoord
    call PRINTLINE(LUVAL_,LINE,124,.false.)
    write(LINE,1010) (nTypes(i),i=1,7)
    call PRINTLINE(LUVAL_,LINE,44,.false.)
    if (isLuscus) then
      write(LINE,1020) iGridNpt(1),iGridNpt(2),iGridNpt(3)
    else
      write(LINE,1020) iGridNpt(1)-1,iGridNpt(2)-1,iGridNpt(3)-1
    end if
    call PRINTLINE(LUVAL_,LINE,20,.false.)
    write(LINE,1030) GridOrigin
    call PRINTLINE(LUVAL_,LINE,47,.false.)
    write(LINE,1040) GridAxis1
    call PRINTLINE(LUVAL_,LINE,47,.false.)
    write(LINE,1050) GridAxis2
    call PRINTLINE(LUVAL_,LINE,47,.false.)
    write(LINE,1060) GridAxis3
    call PRINTLINE(LUVAL_,LINE,47,.false.)
    ! Here we dump all missing information
    if (.not. isLuscus) then
      write(LINE,'(1x,A,I2)') 'CR_SIZE=',iCRSIZE
      call PRINTLINE(LUVAL_,LINE,20,.false.)
      write(LINE,'(1x,A,I2)') 'PACK=',merge(1,0,isMOPack)
      call PRINTLINE(LUVAL_,LINE,20,.false.)
      write(LINE,'(1x,A,I3)') 'BYTES=',NBYTES
      call PRINTLINE(LUVAL_,LINE,20,.false.)
      write(LINE,'(1x,A,I3)') 'N_in_Line=',NINLINE
      call PRINTLINE(LUVAL_,LINE,20,.false.)
      NFIRST = nInc
      if (nBlocks == 1) NFIRST = nCoor
      write(LINE,'(1x,A,I8)') 'N_FIRST=',NFIRST
      call PRINTLINE(LUVAL_,LINE,20,.false.)

      NLAST = nCoor-(nBlocks-1)*nInc
      if (NLAST == 0) NLAST = nInc
      write(LINE,'(1x,A,I8)') 'N_LAST=',NLAST
      call PRINTLINE(LUVAL_,LINE,20,.false.)

      nn = (NFIRST/NINLINE)*(NBYTES*NINLINE+iCRSIZE)
      nn1 = NFIRST-(NFIRST/NINLINE)*NINLINE
      if (nn1 > 0) nn = nn+nn1*NBYTES+iCRSIZE

      write(LINE,'(1x,A,I10)') 'N_OFFSET=',nn
      call PRINTLINE(LUVAL_,LINE,30,.false.)

      nn = (NLAST/NINLINE)*(NBYTES*NINLINE+iCRSIZE)
      nn1 = NLAST-(NLAST/NINLINE)*NINLINE
      if (nn1 > 0) nn = nn+nn1*NBYTES+iCRSIZE
      write(LINE,'(1x,A,I10)') 'N_LAST_OFFSET=',nn
      call PRINTLINE(LUVAL_,LINE,30,.false.)
    end if

    ! skip file pointers
    ! skip end orbital section
  else

    if (isLine) then
      write(line,'(a,2i10)') '# ',nShowMOs_,nCoor
      call PrintLine(LuVal_,line,23,.false.)
    else
      write(line,'(a,a)') 'VERSION=     ',VERSION
      call PrintLine(LuVal_,line,23,.false.)
      !write(line,'(a,a)') 'Extension=   ',0
      !call PrintLine(LuVal_,line,23,.false.)
      write(line,'(a,i10)') 'N_of_MO=     ',nMOs
      !call iXML('nMOs',nMOs)
      call PrintLine(LuVal_,line,23,.false.)
      write(line,'(a,i10)') 'N_of_Grids=  ',nShowMOs_
      !call iXML('nGrids',nShowMOs_)
      call PrintLine(LuVal_,line,23,.false.)
      write(line,'(a,i10)') 'N_of_Points= ',nCoor
      !call iXML('nPoints',nCoor)
      call PrintLine(LuVal_,line,23,.false.)
      write(line,'(a,i10)') 'Block_Size=  ',nInc
      !call iXML('Block Size',nInc)
      call PrintLine(LuVal_,line,23,.false.)
      nBlocks = nCoor/nInc+1
      write(line,'(a,i10)') 'N_Blocks=    ',nBlocks
      !call iXML('nBlocks',nBlocks)
      call PrintLine(LuVal_,line,23,.false.)
      ! new cut off
      write(line,'(a,i10)') 'Is_cutoff=   ',merge(1,0,isCutOff)
      call PrintLine(LuVal_,line,23,.false.)
      write(line,'(a,f10.4)') 'CutOff=      ',CutOff
      call PrintLine(LuVal_,line,23,.false.)
      write(line,'(a,i10)') 'N_P=         ',iiCoord
      call PrintLine(LuVal_,line,23,.false.)
      write(line,'(a,7I5)') 'N_INDEX=     ',nTypes
      call PrintLine(LuVal_,line,48,.false.)
    end if
    if (.not. isTheOne) then

      write(line,'(a,3i5)') 'Net=         ',iGridNpt(1)-1,iGridNpt(2)-1,iGridNpt(3)-1
      !call iaXML('Net',iGridNpt,3)
      call PrintLine(LuVal_,line,28,.false.)

      write(line,'(a,3f12.3)') 'Origin= ',GridOrigin
      !call daXML('Origin',GridOrigin,3)
      call PrintLine(LuVal_,line,44,.false.)

      write(line,'(a,3f12.3)') 'Axis_1= ',GridAxis1
      !call daXML('Axis 1',GridAxis1,3)
      call PrintLine(LuVal_,line,44,.false.)
      write(line,'(a,3f12.3)') 'Axis_2= ',GridAxis2
      !call daXML('Axis 2',GridAxis2,3)
      call PrintLine(LuVal_,line,44,.false.)
      write(line,'(a,3f12.3)') 'Axis_3= ',GridAxis3
      !call daXML('Axis 3',GridAxis3,3)
      call PrintLine(LuVal_,line,44,.false.)

    end if
  end if
end do
!call eXML('Header')

return

1000 format(1X,'N_of_MO=',I5,1X,'N_of_Grids=',I4,1X,'N_of_Points=',I8,1X,'Block_Size=',I6,1X,'N_Blocks=',I4,1X,'Is_cutoff=',I1,1X, &
            'CutOff=',F8.4,1X,'N_P=',I8)
1010 format(1X,'N_INDEX=',7I5)
1020 format(1X,'Net=',3I5)
1030 format(1X,'Origin=',3(1X,F12.8))
1040 format(1X,'Axis_1=',3(1X,F12.3))
1050 format(1X,'Axis_2=',3(1X,F12.3))
1060 format(1X,'Axis_3=',3(1X,F12.3))

end subroutine PrintHeader
