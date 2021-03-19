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
      Subroutine PrintTitles(LuVal,nShowMOs,isDensity,nMOs,             &
     &  iWipGRef, isEner,  WipOcc, iWipType, Crypt,                     &
     &  iWipNZ, WipE, VBocc, ifpartial,isLine,isSphere,                 &
     &  isColor, ISLUSCUS, nCoor,nBlocks,nInc)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
      Implicit Real*8 (A-H,O-Z)

      Character Line*12000
      Character Crypt*7, bb
      Dimension iWipGRef(*), WipOcc(*), iWipType(*),                    &
     &          iWipNZ(*), WipE(*)
       Character LineT*10
       integer Sizeof8
       Sizeof8=8
      iActOrb=0
       LineT='GridName= '
       if(isLine.eq.1) LineT='#GridName='
!       print *,'here',nInc, nBlocks, nCoor
      if(isLUSCUS.eq.1) then
!        NFIRST=nInc

!        if(nBlocks.eq.1) NFIRST=nCoor

!      ii=0
!      do i=1,nShowMOs
!      write(Line,'(A17,1000i12)') ' File_pointers = ',
!     * ((nFirst*(j-1)+ii)*Sizeof8,j=1,nBlocks)
!      CALL PRINTLINE(LUVAL,LINE,12*nBlocks+17,0)
!      ii=ii+nFirst*nBlocks
!      enddo
      write(Line,'(A,i22)') ' ORBOFF = ', nCoor*nShowMOs*Sizeof8
!     * + 28*nShowMOs - 20
! ?? -20?
      CALL PRINTLINE(LUVAL,LINE,32,0)
      endif
      do i=1,nShowMOs-isDensity-isSphere-isColor
        j=iWipGRef(i)
        if(isEner.eq.1) then
          if(.not.(0.eq.1.and.                                          &
     &       WipOcc(j).gt.0d0.and.                                      &
     &       WipOcc(j).lt.2d0)) then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            IF (ISLUSCUS .EQ. 1) THEN
              WRITE(LINE,1000) iWipNZ(j), iWipNZ(j+nMOs), WipE(j),      &
     &                         WipOcc(j), bb
              CALL PRINTLINE(LUVAL,LINE,72,0)
            ELSE
              write (line,'(a,i2,i5,f12.4,'' ('',f4.2,'')'',1x,a)')     &
     &                    LineT,                                        &
     &                    iWipNZ(j),iWipNZ(j+nMOs),                     &
     &                    WipE(j),WipOcc(j),bb
              call PrintLine(LuVal,line,38,0)
            END IF
 1000       FORMAT(1X,'GridName= Orbital sym=',i2,' index=',i5,         &
     &             ' Energ=',F12.4,' occ=', F4.2,' type=',a1)
          else
            iActOrb=iActOrb+1
            IF (ISLUSCUS .EQ. 1) THEN
              WRITE(LINE,1010) 'VB orbital',iActOrb,' (',VBocc,')'
 1010         FORMAT(1X,'GridName= VB_orbital iActOrb= ',I4,            &
     &               ' occ= ',F4.2)
              CALL PRINTLINE(LUVAL,LINE,45,0)
            ELSE
              write (line,'(2a,i4,5x,a,f4.2,a)') LineT,                 &
     &                    'VB orbital',iActOrb,' (',VBocc,')'
              call PrintLine(LuVal,line,38,0)
            END IF
          endif
        else
          if(.not.(0.eq.1.and.                                          &
     &       WipOcc(j).gt.0d0 .and.                                     &
     &       WipOcc(j).lt.2d0)) then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            IF (ISLUSCUS .EQ. 1) THEN
              WRITE(LINE,1020) iWipNZ(j), iWipNZ(j+nMOs),               &
     &                         WipOcc(j), bb
 1020         FORMAT(1X,'GridName= Orbital sym=',I2,' index=',I5,       &
     &                  ' occ=',F4.2,' type=',A1)
              CALL PRINTLINE(LUVAL,LINE,53,0)
            ELSE
              write (line,'(a,i2,i5,'' ('',f8.6,'')'',1x,a)')           &
     &                    LineT,iWipNZ(j),                              &
     &                    iWipNZ(j+nMOs), WipOcc(j),bb
              call PrintLine(LuVal,line,30,0)
            END IF
          else
            iActOrb=iActOrb+1
            IF (ISLUSCUS .EQ. 1) THEN
              WRITE(LINE,1010) iActOrb, VBocc
              CALL PRINTLINE(LUVAL,LINE,45,0)
            ELSE
              write (line,'(2a,i4,5x,a,f4.2,a)') LineT,                 &
     &                         'VB orbital',iActOrb,' (',VBocc,')'
              call PrintLine(LuVal,line,30,0)
            END IF
          endif
        endif
      enddo
      if(isSphere.eq.1) Then
          write (line,'(a,a)') LineT,'  Sphere '
          call PrintLine(LuVal,line,19,0)
      endif
      if(isSphere.eq.1) Then
          write (line,'(a,a)') LineT,'0 Color  '
          call PrintLine(LuVal,line,19,0)
      endif
      if(isDensity.eq.1) Then
        if(ifpartial.eq.0) Then
          IF (ISLUSCUS .EQ. 1) THEN
            LINE=' GridName= Density'
            CALL PRINTLINE(LUVAL,LINE,18,0)
          ELSE
            write (line,'(a,a)') LineT,'  Density'
            call PrintLine(LuVal,line,19,0)
          END IF
        else
          IF (ISLUSCUS .EQ. 1) THEN
            LINE=' GridName= Density (partial)'
            CALL PRINTLINE(LUVAL,LINE,28,0)
          ELSE
            write (line,'(a,a)') LineT,'  Density (partial)'
            call PrintLine(LuVal,line,29,0)
          END IF
        endif
      endif
          IF (ISLUSCUS .EQ. 1) THEN
            LINE=' <DENSITY>'
            CALL PRINTLINE(LUVAL,LINE,10,0)
          endif

      return
! Avoid unused argumet warnings
      if (.false.) then
        call Unused_integer(nBlocks)
        call Unused_integer(nInc)
      end if
      end
