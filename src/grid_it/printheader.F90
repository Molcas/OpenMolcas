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
      SubRoutine PrintHeader(nMOs,nShowMOs,nShowMOs_ab,nCoor,nInc,      &
     &  iiCoord,nTypes,iCRSIZE,NBYTES,NINLINE,nBlocks )
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "grid.fh"
      character line*128
      Integer nTypes(7)
      LuVal_=LuVal
      if(isLUSCUS.eq.1) LuVal_=LID
      nShowMOs_=nShowMOs
!      call bXML('Header')
!      call iXML('UHF',isUHF)
      do iiUHF=0,isUHF
      if(iiUHF.eq.1) then
      LuVal_=LuVal_ab
      if(isLUSCUS.eq.1) LuVal_=LID_ab
      nShowMOs_=nShowMOs_ab
      endif
      IF (ISLUSCUS .EQ. 1) THEN
! debug
!        write(*,*) 'N_of_MO=', nMOs
!        write(*,*) 'N_of_Grids=', nShowMOs_
!        write(*,*) 'N_of_Points=', nCoor
!        write(*,*) 'Block_Size=', nInc
!        write(*,*) 'N_Blocks=', nBlocks
!        write(*,*) 'Is_cutoff=', isCutOff
!        write(*,*) 'CutOff=', Cutoff
!        write(*,*) 'N_P=', iiCoord
! end debug


        IF (nMOs .GT. 99999) THEN
          Write(6,*) 'Number of MO''s can''t be larger that 99999'
          Call Quit_OnUserError()
        END IF
        IF (nShowMOs_ .GT. 9999) THEN
          Write(6,*) 'Number of grids can''t be larger that 9999'
          Call Quit_OnUserError()
        END IF
        IF (nCoor .GT. 99999999) THEN
          Write(6,*) 'Number of points can''t be larger that 99999999'
          Call Quit_OnUserError()
        END IF
        IF (nInc .GT. 999999) THEN
          Write(6,*) 'Block size can''t be larger that 999999'
          Call Quit_OnUserError()
        END IF
        nBlocks=nCoor/nInc+1
        IF (nBlocks .GT. 9999) THEN
          Write(6,*) 'Number of blocks can''t be larger that 9999'
!          write(*,*) 'nBlocks = ', nBlocks
          Call Quit_OnUserError()
        END IF
        IF (isCutOff .GT. 9) THEN
          Write(6,*) 'Wrong cutoff'
          Call Quit_OnUserError()
        END IF
        IF (iiCoord .GT. 99999999) THEN
          Write(6,*) 'N_P can''t be larger that 99999999'
          Call Quit_OnUserError()
        END IF
        WRITE(LINE,'(''<GRID>'')')
        CALL PRINTLINE(LUVAL_, LINE,6,0)

        WRITE(LINE,1000)                                                &
     &         nMOs, nShowMOs_, nCoor, nInc, nBlocks,                   &
     &         isCutOff, Cutoff, iiCoord
 1000   FORMAT(1X,'N_of_MO=',I5,1X,'N_of_Grids=',I4,                    &
     &         1X,'N_of_Points=',I8,1X,'Block_Size=',I6,                &
     &         1X,'N_Blocks=',I4,1X,'Is_cutoff=',I1,                    &
     &         1X,'CutOff=',F8.4,1X,'N_P=',I8)
        CALL PRINTLINE(LUVAL_, LINE, 124,0)
        WRITE(LINE,1010) (nTypes(i),i=1,7)
 1010   FORMAT(1X,'N_INDEX=',7I5)
        CALL PRINTLINE(LUVAL_, LINE, 44,0)
        if(isLUSCUS.eq.0) then
        WRITE(LINE,1020) iGridNpt(1)-1, iGridNpt(2)-1, iGridNpt(3)-1
        else
        WRITE(LINE,1020) iGridNpt(1), iGridNpt(2), iGridNpt(3)
        endif
 1020   FORMAT(1X,'Net=',3I5)
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        WRITE(LINE,1030) GridOrigin
 1030   FORMAT(1X,'Origin=',3(1X,F12.8))
        CALL PRINTLINE(LUVAL_, LINE, 47,0)
        WRITE(LINE,1040) GridAxis1
 1040   FORMAT(1X,'Axis_1=',3(1X,F12.3))
        CALL PRINTLINE(LUVAL_, LINE, 47,0)
        WRITE(LINE,1050) GridAxis2
 1050   FORMAT(1X,'Axis_2=',3(1X,F12.3))
        CALL PRINTLINE(LUVAL_, LINE, 47,0)
        WRITE(LINE,1060) GridAxis3
 1060   FORMAT(1X,'Axis_3=',3(1X,F12.3))
        CALL PRINTLINE(LUVAL_, LINE, 47,0)
! Here we dump all missing information
        if(isLUSCUS.eq.0) then
        WRITE(LINE,'(1x,A,I2)') 'CR_SIZE=',iCRSIZE
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        WRITE(LINE,'(1x,A,I2)') 'PACK=',imoPACK
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        WRITE(LINE,'(1x,A,I3)') 'BYTES=',NBYTES
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        WRITE(LINE,'(1x,A,I3)') 'N_in_Line=',NINLINE
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        NFIRST=nInc
        if(nBlocks.eq.1) NFIRST=nCoor
        WRITE(LINE,'(1x,A,I8)') 'N_FIRST=',NFIRST
        CALL PRINTLINE(LUVAL_, LINE, 20,0)

        NLAST=nCoor-(nBlocks-1)*nInc
        if(NLAST.eq.0) NLAST=nInc
        WRITE(LINE,'(1x,A,I8)') 'N_LAST=',NLAST
        CALL PRINTLINE(LUVAL_, LINE, 20,0)

        nn=(NFIRST/NINLINE)*(NBYTES*NINLINE+iCRSIZE)
        nn1=NFIRST-(NFIRST/NINLINE)*NINLINE
        if(nn1.gt.0) nn=nn+nn1*NBYTES+iCRSIZE

        WRITE(LINE,'(1x,A,I10)') 'N_OFFSET=',nn
        CALL PRINTLINE(LUVAL_, LINE, 30,0)

        nn=(NLAST/NINLINE)*(NBYTES*NINLINE+iCRSIZE)
        nn1=NLAST-(NLAST/NINLINE)*NINLINE
        if(nn1.gt.0) nn=nn+nn1*NBYTES+iCRSIZE
        WRITE(LINE,'(1x,A,I10)') 'N_LAST_OFFSET=',nn
        CALL PRINTLINE(LUVAL_, LINE, 30,0)
        endif

!
! skip file pointers
! skip end orbital section
!
      ELSE


      if(isLine.eq.0) then
       write (line,'(a,a)') 'VERSION=     ',VERSION
       call PrintLine(LuVal_,line,23,0)
!       write (line,'(a,a)') 'Extension=   ',0
!       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)') 'N_of_MO=     ',nMOs
!       call iXML('nMOs',nMOs)
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)') 'N_of_Grids=  ',nShowMOs_
!       call iXML('nGrids',nShowMOs_)
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)') 'N_of_Points= ',nCoor
!       call iXML('nPoints',nCoor)
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)') 'Block_Size=  ',nInc
!       call iXML('Block Size',nInc)
       call PrintLine(LuVal_,line,23,0)
       nBlocks=nCoor/nInc+1
       write (line,'(a,i10)')   'N_Blocks=    ',nBlocks
!       call iXML('nBlocks',nBlocks)
       call PrintLine(LuVal_,line,23,0)
! new cut off
       write (line,'(a,i10)')   'Is_cutoff=   ',isCutOff
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,f10.4)') 'CutOff=      ',CutOff
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)')   'N_P=         ',iiCoord
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,7I5)')   'N_INDEX=     ',nTypes
       call PrintLine(LuVal_,line,48,0)
      else
       write(line,'(a,2i10)') '# ', nShowMOs_ , nCoor
       call PrintLine(LuVal_,line,23,0)
      endif
      if (isTheOne .eq. 1) goto 777

      write (line,'(a,3i5)') 'Net=         ',iGridNpt(1)-1,             &
     &        iGridNpt(2)-1,iGridNpt(3)-1
!      call iaXML('Net',iGridNpt,3)
      call PrintLine(LuVal_,line,28,0)

      write (line,'(a,3f12.3)') 'Origin= ',GridOrigin
!      call daXML('Origin',GridOrigin,3)
      call PrintLine(LuVal_,line,44,0)

      write (line,'(a,3f12.3)') 'Axis_1= ',GridAxis1
!      call daXML('Axis 1',GridAxis1,3)
      call PrintLine(LuVal_,line,44,0)
      write (line,'(a,3f12.3)') 'Axis_2= ',GridAxis2
!      call daXML('Axis 2',GridAxis2,3)
      call PrintLine(LuVal_,line,44,0)
      write (line,'(a,3f12.3)') 'Axis_3= ',GridAxis3
!      call daXML('Axis 3',GridAxis3,3)
      call PrintLine(LuVal_,line,44,0)

777   continue
      END IF
      enddo
!      call eXML('Header')

      Return
      End
