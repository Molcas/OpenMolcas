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
      Subroutine DumpM2Msi(iRun,Luval,LID,nShowMOs,isDensity,nMOs,      &
     &  iWipGRef,  WipOcc, WipMO, WipOut, mCoor,                        &
     &   iGauss, nInc, imoPack, iWipPBlock,                             &
     &  cMoBlock,nBytesPackedVal, dnorm, Crypt, VbOcc,                  &
     &  isTheOne,isLine,isBinary, isEner, iWipType, iWipNZ,WipE,        &
     &  WLine,nLine,WCoor,iPrintCount,isDebug,                          &
     &  isCutOff, iWipCutOff,isSphere,SphrDist,isColor,SphrColor,       &
     &  ISLUSCUS, NBYTES,NINLINE)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"
!      Dimension xLimits(4)
!      Integer iYDelta(3)
      Character*1  cMoBlock(*)
      Character Crypt*7, bb
      Character Line*128
!     Character fmt*20
! test
!     character*3 cint
!     character*1 cx(64)
!     data cx /
!    +         '0','1','2','3','4','5','6','7','8','9',
!    +         'a','b','c','d','e','f','g','h','i','j',
!    +         'k','l','m','n','o','p','q','r','s','t',
!    +         'u','v','w','x','y','z','A','B','C','D',
!    +         'E','F','G','H','I','J','K','L','M','N',
!    +         'O','P','Q','R','S','T','U','V','W','X',
!    +         'Y','Z','@','#' /
! test end
      Dimension iWipGRef(*), WipOcc(*), WipMO(*), WipOut(*),            &
     &   iWipPBlock(*),iWipType(*), iWipNZ(*),WipE(*),                  &
     &  WLine(nLine,mCoor),WCoor(3,mCoor),iWipCutOff(*),                &
     &  SphrDist(mCoor),SphrColor(mCoor)
       Dimension DumArr(2)
#include "WrkSpc.fh"
!          write(*,*) 'entering DumpM2Msi'
           if(irun.gt.100) print *, iGauss, nbytes, ninc, ninline
          iActOrb=0
          iPrintCount=iPrintCount+1
        do i=1, nShowMOs-isDensity-isSphere-isColor
          iMOs=iWipGRef(i)

          if(.not.(0.eq.1.and.                                          &
     &       WipOcc(iMOs).gt.0d0                                        &
     &        .and.                                                     &
     &       WipOcc(iMOs).lt.2d0)) then
            call outmo(iMOs,1,WipMO,dumArr,WipOut,mCoor,nMOs)
!          else
!            iActOrb=iActOrb+1
!            call outmo(0,1,WipMO,WipVBmat(1+(iActOrb-1)*nMOs),
!     >           WipOut,mCoor,nMOs)
          endif

        if(isLine.eq.0.and.isLuscus.eq.0) then
            write (line,'(a,i4)') 'Title= ',iMOs
            call PrintLine(LuVal,line,12,1)
          endif
          if(isTheOne.eq.1)then
            if(isLine.eq.0 .AND. ISLUSCUS .eq. 0) then
              write(LuVal,'(f18.12)')                                   &
     &            (WipOut(j),j=1,mCoor)
            else
              if(i+1.le.nLine) then
                do j=1,mCoor
                  WLine(i+1,j)=WipOut(j)
                enddo
              endif
            endif
            goto 3939
          endif

!          if (imoPack .ne. 0) then
!c            Call PackBlock(WipOut,iWipPBlock,mCoor,
!c     >                     xLimits,iYDelta)
!            write (line,9000)
!     *                     0,
!     *                     (xLimits(j),j=1,4),(iYDelta(j),j=1,3)
!            call PrintLine(LuVal,line,73,0)
!
!9000        format ('BHeader=',I2,1X,(4(E10.4,1X),3(I5,1X)))
!            if (isBinary .ne. 0) then
!c              call IArrToChar(iWipPBlock,cMoBlock,mCoor)
!cvv ! NOT CODED YET
!              write (LuVal) (cMoBlock(j),j=1,mCoor*nBytesPackedVal)
!            else
!              write (LuVal,'(I5)') (iWipPBlock(j), j=0,mCoor-1)
!            endif
!          else
            if (isBinary .eq. 1) then
! packing late
              if(isCutOff.eq.1) then
                Call GetMem('TMP','ALLO','REAL',ipCMP,mCoor)
                call dcopy_(mCoor,WipOut,1,Work(ipCMP),1)
                iii=0
                do ii=1,mCoor
                  if(iWipCutOff(ii).eq.1) then
                   Work(ipCMP+iii)=WipOut(ii)
                   iii=iii+1
                  endif
                enddo

                write (LuVal) (Work(ipCMP+j),j=0,iii-1)
                Call GetMem('TMP','FREE','REAL',ipCMP,mCoor)
              else
! no cut off
                write (LuVal) (WipOut(j),j=1,mCoor)

              endif
            else !isBinary
!             write(cint, '(i3.3)') i
              if(isDebug.eq.0) then
! normal output - just numbers
               if(isCutOff.eq.1) then
                  do j=1,mCoor
                    if(iWipCutOff(j).eq.1)                              &
     &                write (LuVal, '(E10.4)') WipOut(j)
                  enddo
               else if (isLUSCUS .eq. 1) then
! NOPACKING
!                if(imoPack.eq.0) then
                 call dump_lusc(LID, WipOut,mCoor)
! debug dump of data
!
!                  do j=1,mCoor,NINLINE
!                    num=mCoor-j
!                    if(num.gt.NINLINE) num=NINLINE
!                    write(fmt,'(A,I2,A,I2,A,A)')
!     &                       '(',NINLINE,'E',NBYTES,'.4',')'
!                    write(line, fmt) (WipOut(j-1+ij),ij=1,num)
!                    call printline(LID,line,NINLINE*NBYTES,0)
!                  enddo
!                else
!c PACKING NOT implemented
!C                  open(unit=38, file='testgr_'//cint//'.txt')
!                   do j = 1, min(mcoor,10)
!                     iexpnt=int(log10(abs(wipout(j))))
!                     write(6,*) 'num=', wipout(j), 'exp = ', iexpnt
!                     dnum=(1.0D0 + wipout(j)/10.0D0**iexpnt)/2.0D0
!
!                     in1 = int(dnum*64.0d0)
!                     in2 = int((dnum-dble(in1)*64.0D0)*4096.0D0)
!                     iexpnt = iexpnt + 50
!                     write(6,'(1x,3(1x,e18.8),2x,3(1x,i3))')
!     +                      wipout(j), dnum, dexpnt,
!     +                      in1, in2, iexpnt
!                     if (iexpnt .lt. 1) then
!                       iexpnt=1
!                     else if (iexpnt .gt. 64) then
!                       iexpnt=64
!                     end if
!                     write(6,'(1x,3(1x,e18.8),2x,3(1x,i3))')
!     +                      wipout(j), dnum, dexpnt,
!     +                      in1, in2, iexpnt
!c                     write(*,*) '-----------------------'
!
!C     +                      cx(in1), cx(in2), cx(iexpnt)
!                   end do
!C                  close(38)
!c                  xxxmin=9.99d+99
!c                  xxxmax=-9.99d+99
!c                  do j=1,mCoor
!c                    if(wipout(j) .gt. xxxmax) xxxmax = wipout(j)
!c                    if(wipout(j) .lt. xxxmin) xxxmin = wipout(j)
!c                  end do
!C                  write(*, *) 'test min/max', xxxmin, xxxmax
!
!
!                endif !imoPack
               else!isCutOff
! writing of data
                  write (LuVal,'(E10.4)')                               &
     &                      (WipOut(j),j=1,mCoor)
               endif !isCutOff, isLuscus
              else !isDebug
! extended output -
                write (LuVal,'(E10.4,3f8.4)')                           &
     &       (WipOut(j),WCoor(1,j), WCoor(2,j), Wcoor(3,j) ,j=1,mCoor)
              endif !isDebug
            endif !isBinary
!          endif !imoPack

        j=iWipGRef(i)

        if(isEner.eq.1) then
          if(.not.(0.eq.1.and.                                          &
     &       WipOcc(j).gt.0d0.and.WipOcc(j).lt.2d0))then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            if(iRun.eq.1.and.iPrintCount.eq.1)                          &
     &      write(6,'(a,i2,i5,f12.4,'' ('',f4.2,'') '',a)')             &
     &            'GridName= ',                                         &
     &            iWipNZ(j),iWipNZ(j+nMOs),                             &
     &            WipE(j),WipOcc(j),bb
          else
!            iActOrb=iActOrb+1
            if(iRun.eq.1.and.iPrintCount.eq.1)                          &
     &      write(6,'(2a,i4,5x,a,f4.2,a)') 'GridName= ',                &
     &            'VB orbital',iActOrb,' (',VBocc,')'
          endif
        else
          if(.not.(0.eq.1.and.                                          &
     &       WipOcc(j).gt.0d0.and.WipOcc(j).lt.2d0))then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            if(iRun.eq.1.and.iPrintCount.eq.1)                          &
     &      write(6,'(a,i2,i5,'' ('',f8.6,'') '',a)')                   &
     &            'GridName= ',                                         &
     &            iWipNZ(j),iWipNZ(j+nMOs),                             &
     &            WipOcc(j),bb
          else
!            iActOrb=iActOrb+1
            if(iRun.eq.1.and.iPrintCount.eq.1)                          &
     &      write(6,'(2a,i4,5x,a,f4.2,a)') 'GridName= ',                &
     &            'VB orbital',iActOrb,' (',VBocc,')'
          endif
        endif
3939      continue
        enddo

        if(isSphere.eq.1) then
!VV FIXME: no packing.
            write (line,'(a,i4)') 'Title= ',-1
            call PrintLine(LuVal,line,12,1)
          if(isBinary.eq.0) then
           do j=1,mCoor
              write (LuVal,'(E18.12)')                                  &
     &        SphrDist(j)
           enddo
          else
              write(LuVal) (SphrDist(j),j=1,mCoor)
          endif
        endif

        if(isColor.eq.1) then
!VV FIXME: no packing.
            write (line,'(a,i4)') 'Title= ',-10
            call PrintLine(LuVal,line,12,1)
          if(isBinary.eq.0) then
           do j=1,mCoor
              write (LuVal,'(E18.12)')                                  &
     &        SphrColor(j)
           enddo
          else
              write(LuVal) (SphrColor(j),j=1,mCoor)
          endif
        endif

        if(isDensity.eq.1) then
          call outmo(0,2,WipMO,WipOcc,WipOut,                           &
     &               mCoor,nMOs)
          do j=1, mCoor
            dNorm=dNorm+WipOut(j)
          enddo
!****
!          write(6,*) " mCoor=",mCoor
!****
          if(isLine.eq.0.and.IsLuscus.eq.0) then
            write (line,'(a,i4)') 'Title= ',0
            call PrintLine(LuVal,line,12,1)
          endif
!          if (imoPack.ne.0) then
!c         print *,'pack code'
!c            Call PackBlock(WipOut,iWipPBlock,mCoor,
!c     >                     xLimits,iYDelta)
!c            write (line,9000)
!c     *                    0,
!c     *                    (xLimits(j),j=1,4),(iYDelta(j),j=1,3)
!c            call PrintLine(LuVal,line,73,0)
!c
!c            if (isBinary .ne. 0) then
!c              call IArrToChar(iWipPBlock,cMoBlock,mCoor)
!cvv!!!!!!!
!c              IF (ISLUSCUS .EQ. 1) THEN
!c                RC=C_WRITE(LID, CMOBLOCK, (mCoor*nBytesPackedVal)*RTOB) !!!!!!!!!!!!!!!!!!!!check mCoor*nBytesPackedVal
!c                IF (RC .EQ. 0) THEN
!c                  WRITE(6,*) 'error in writing luscus file!'
!c                  CALL Abend()
!c                END IF
!c              ELSE
!c                write (LuVal) (cMoBlock(j),j=1,mCoor*nBytesPackedVal)
!c              END IF
!c            else
!c              IF (ISLUSCUS .EQ. 1) THEN
!c                RC=C_WRITE(LID, IWIPPBLOCK, (mCoor)*RTOB) !!!!!!!!!!!!!!!!!!!!check mCoor*nBytesPackedVal
!c                IF (RC .EQ. 0) THEN
!c                  WRITE(6,*) 'error in writing luscus file!'
!c                  CALL Abend()
!c                END IF
!c              ELSE
!c                write (LuVal,'(I5)') (iWipPBlock(j), j=1,mCoor)
!c              END IF
!c            endif
!          else
              IF (ISLUSCUS .EQ. 1) THEN
                 call dump_lusc(LID, WipOut,mCoor)
               Endif

            if (isBinary .eq. 1) then
! packing late
              if(isCutOff.eq.1) then
                Call GetMem('TMP','ALLO','REAL',ipCMP,mCoor)
                call dcopy_(mCoor,WipOut,1,Work(ipCMP),1)
                iii=0
                do ii=1,mCoor
                  if(iWipCutOff(ii).eq.1) then
                   Work(ipCMP+iii)=WipOut(ii)
                   iii=iii+1
                  endif
                enddo

            IF (ISLUSCUS .EQ. 1) THEN
              !!!!!!!!!!!!!!!!!!!!check iii-1
              RC=C_WRITE(LID, WORK(IPCMP), (III-1)*RTOB)
              IF (RC .EQ. 0) THEN
                WRITE(6,*) 'error in writing luscus file!'
                CALL Abend()
              END IF
            ELSE
              write (LuVal) (Work(ipCMP+j),j=0,iii-1)
            END IF
            Call GetMem('TMP','FREE','REAL',ipCMP,mCoor)
           else
! no cut off
!              IF (ISLUSCUS .EQ. 1) THEN
!                 call dump_lusc(LID, WipOut,mCoor)
!           print *,'here'
!                RC=C_WRITE(LID, WIPOUT, MCOOR*RTOB) !!!!!!!!!!!!!!!!!!!!check MCOOR
!                IF (RC .EQ. 0) THEN
!                  WRITE(6,*) 'error in writing luscus file!'
!                  CALL Abend()
!                END IF
!              ELSE
                 write (LuVal) (WipOut(j),j=1,mCoor)
!              END IF
            endif
            else

               if(isLine.eq.1) then
                do j=1,mCoor
                  WLine(1,j)=WipOut(j)
                enddo
               else
                 if (isDebug.eq.0) then
! normal output - just numbers
            if(isCutOff.eq.1) then
              do j=1,mCoor
                if(iWipCutOff(j).eq.1)                                  &
     &            write (LuVal, '(E18.12)') WipOut(j)
              enddo
             else

              if(ISLUSCUS.eq.0)                                         &
     &              write (LuVal,'(E18.12)') (WipOut(j),j=1,mCoor)
              endif
                 else
! extra output -
              if(ISLUSCUS.eq.0) write (LuVal,'(E18.12,3f8.4)')          &
     &       (WipOut(j),WCoor(1,j), WCoor(2,j), Wcoor(3,j) ,j=1,mCoor)
                 endif

               endif
!GG This is only for testing CASDFT functional. It will be restore.
!GG              write (LuVal,'(E10.4)') (Work(j),j=ipOut,ipOut+mCoor-1)
            endif
          endif
!        endif

      if(isLine.eq.1.and.ISLUSCUS.eq.0) then
       do i=1,mCoor
        write(LuVal,'(3F10.6,22E20.12)')                                &
     &  (WCoor(j,i),j=1,3),(WLine(j,i),j=1,nLine)
       enddo
      endif

!
      Return
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(imoPack)
         Call Unused_integer_array(iWipPBlock)
         Call Unused_character_array(cMoBlock)
         Call Unused_integer(nBytesPackedVal)
      End If
      end
