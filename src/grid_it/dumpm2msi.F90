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

subroutine DumpM2Msi(iRun,Luval,LID,nShowMOs,isDensity,nMOs,iWipGRef,WipOcc,WipMO,WipOut,mCoor,iGauss,nInc,isMOPack,iWipPBlock, &
                     cMoBlock,nBytesPackedVal,dnorm,Crypt,VbOcc,isTheOne,isLine,iBinary,isEner,iWipType,iWipNZ,WipE,WLine,nLine, &
                     WCoor,iPrintCount,isDebug,isCutOff,iWipCutOff,isSphere,SphrDist,isColor,SphrColor,isLuscus,NBYTES,NINLINE)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Ten
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: iRun, LuVal, LID, nShowMOs, nMOs, iWipGRef(*), mCoor, iGauss, nInc, iWipPBlock(*), &
                                 nBytesPackedVal, iBinary, iWipType(*), iWipNZ(*), nLine, iWipCutOff(*), NBYTES, NINLINE
logical(kind=iwp), intent(in) :: isDensity, isMOPack, isTheOne, isLine, isEner, isDebug, isCutOff, isSphere, isColor, isLuscus
integer(kind=iwp), intent(inout) :: iPrintCount
real(kind=wp), intent(in) :: WipOcc(*), WipMO(*), WipOut(*), VbOcc, WipE(*), WCoor(3,mCoor), SphrDist(mCoor), SphrColor(mCoor)
character, intent(in) :: cMoBlock(*)
real(kind=wp), intent(inout) :: dNorm, WLine(nLine,mCoor)
character(len=7), intent(in) :: Crypt
integer(kind=iwp) :: i, iActOrb, ib, ii, iii, iMOs, j, RC !, iYDelta(3)
real(kind=wp) :: DumArr(2) !, xLimits(4)
character :: bb
character(len=128) :: Line
real(kind=wp), allocatable :: CMP(:)
!character(len=20) :: formt
! test
!character(len=3) :: cint
!character, parameter :: cx(64) = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ@#'
! test end
integer(kind=iwp), external :: C_WRITE

!write(u6,*) 'entering DumpM2Msi'
if (irun > 100) print *, iGauss,nbytes,ninc,ninline
iActOrb = 0
iPrintCount = iPrintCount+1
do i=1,nShowMOs-merge(1,0,isDensity)-merge(1,0,isSphere)-merge(1,0,isColor)
  iMOs = iWipGRef(i)

  if (.not.(.false. .and. (WipOcc(iMOs) > Zero) .and. (WipOcc(iMOs) < Two))) then
    call outmo(iMOs,1,WipMO,dumArr,WipOut,mCoor,nMOs)
  !else
  !  iActOrb = iActOrb+1
  !  call outmo(0,1,WipMO,WipVBmat(1+(iActOrb-1)*nMOs),WipOut,mCoor,nMOs)
  end if

  if ((.not. isLine) .and. (.not. isLuscus)) then
    write(line,'(a,i4)') 'Title= ',iMOs
    call PrintLine(LuVal,line,12,.true.)
  end if
  if (isTheOne) then
    if ((.not. isLine) .and. (.not. isLuscus)) then
      write(LuVal,'(f18.12)') (WipOut(j),j=1,mCoor)
    else
      if (i+1 <= nLine) then
        do j=1,mCoor
          WLine(i+1,j) = WipOut(j)
        end do
      end if
    end if
  else

    !if (isMOPack) then
    !  !call PackBlock(WipOut,iWipPBlock,mCoor,xLimits,iYDelta)
    !  write(line,9000) 0,(xLimits(j),j=1,4),(iYDelta(j),j=1,3)
    !  call PrintLine(LuVal,line,73,.false.)
    !  9000 format ('BHeader=',I2,1X,(4(E10.4,1X),3(I5,1X)))
    !  if (iBinary /= 0) then
    !    !call IArrToChar(iWipPBlock,cMoBlock,mCoor)
    !    !vv ! NOT CODED YET
    !    write(LuVal) (cMoBlock(j),j=1,mCoor*nBytesPackedVal)
    !  else
    !    write(LuVal,'(I5)') (iWipPBlock(j), j=0,mCoor-1)
    !  end if
    !else
    if (iBinary == 1) then
      ! packing late
      if (isCutOff) then
        call mma_allocate(CMP,mCoor,label='TMP')
        iii = 0
        do ii=1,mCoor
          if (iWipCutOff(ii) == 1) then
            iii = iii+1
            CMP(iii) = WipOut(ii)
          end if
        end do

        write(LuVal) CMP(1:iii)
        call mma_deallocate(CMP)
      else
        ! no cut off
        write(LuVal) (WipOut(j),j=1,mCoor)

      end if
    else !iBinary
      !write(cint,'(i3.3)') i
      if (isDebug) then
        ! extended output -
        write(LuVal,'(E10.4,3f8.4)') (WipOut(j),WCoor(1,j),WCoor(2,j),Wcoor(3,j),j=1,mCoor)
      else !isDebug
        ! normal output - just numbers
        if (isCutOff) then
          do j=1,mCoor
            if (iWipCutOff(j) == 1) write(LuVal,'(E10.4)') WipOut(j)
          end do
        else if (isLUSCUS) then
          ! NOPACKING
          !if (.not. isMOPack) then
          call dump_lusc(LID,WipOut,mCoor)
          ! debug dump of data
          !
          !  do j=1,mCoor,NINLINE
          !    num = mCoor-j
          !    if (num > NINLINE) num = NINLINE
          !    write(formt,'(A,I2,A,I2,A,A)') '(',NINLINE,'E',NBYTES,'.4',')'
          !    write(line,formt) (WipOut(j-1+ij),ij=1,num)
          !    call printline(LID,line,NINLINE*NBYTES,.false.)
          !  end do
          !else
          !  ! PACKING NOT implemented
          !  !open(unit=38,file='testgr_'//cint//'.txt')
          !  do j=1,min(mcoor,10)
          !    iexpnt = int(log10(abs(WipOut(j))))
          !    write(u6,*) 'num=',WipOut(j),'exp = ',iexpnt
          !    dnum = (One+WipOut(j)/Ten**iexpnt)/Two
          !
          !    in1 = int(dnum*64.0_wp)
          !    in2 = int((dnum-in1*64.0_wp)*4096.0_wp)
          !    iexpnt = iexpnt+50
          !    write(u6,'(1x,3(1x,e18.8),2x,3(1x,i3))') WipOut(j),dnum,dexpnt,in1,in2,iexpnt
          !    if (iexpnt < 1) then
          !      iexpnt = 1
          !    else if (iexpnt > 64) then
          !      iexpnt = 64
          !    end if
          !    write(u6,'(1x,3(1x,e18.8),2x,3(1x,i3))') WipOut(j),dnum,dexpnt,in1,in2,iexpnt
          !    !                                       cx(in1:in1),cx(in2:in2),cx(iexpnt:iexpnt)
          !    !write(u6,*) '-----------------------'
          !  end do
          !  !close(38)
          !  !xxxmin = huge(xxxmin)
          !  !xxxmax = -huge(xxxmax)
          !  !do j=1,mCoor
          !  !  if (WipOut(j) > xxxmax) xxxmax = WipOut(j)
          !  !  if (WipOut(j) < xxxmin) xxxmin = WipOut(j)
          !  !end do
          !  !write(u6,*) 'test min/max',xxxmin,xxxmax
          !
          !end if !isMOPack
        else !isCutOff
          ! writing of data
          write(LuVal,'(E10.4)') (WipOut(j),j=1,mCoor)
        end if !isCutOff, isLuscus
      end if !isDebug
    end if !iBinary
    !end if !isMOPack

    j = iWipGRef(i)

    if (isEner) then
      if (.not.(.false. .and. (WipOcc(j) > Zero) .and. (WipOcc(j) < Two))) then
        ib = iWipType(j)
        bb = ' '
        if ((ib > 0) .and. (ib < 8)) bb = Crypt(ib:ib)
        if ((iRun == 1) .and. (iPrintCount == 1)) &
          write(u6,'(a,i2,i5,f12.4," (",f4.2,") ",a)') 'GridName= ',iWipNZ(j),iWipNZ(j+nMOs),WipE(j),WipOcc(j),bb
      else
        !iActOrb = iActOrb+1
        if ((iRun == 1) .and. (iPrintCount == 1)) write(u6,'(2a,i4,5x,a,f4.2,a)') 'GridName= ','VB orbital',iActOrb,' (',VBocc,')'
      end if
    else
      if (.not.(.false. .and. (WipOcc(j) > Zero) .and. (WipOcc(j) < Two))) then
        ib = iWipType(j)
        bb = ' '
        if ((ib > 0) .and. (ib < 8)) bb = Crypt(ib:ib)
        if ((iRun == 1) .and. (iPrintCount == 1)) &
          write(u6,'(a,i2,i5," (",f8.6,") ",a)') 'GridName= ',iWipNZ(j),iWipNZ(j+nMOs),WipOcc(j),bb
      else
        !iActOrb = iActOrb+1
        if ((iRun == 1) .and. (iPrintCount == 1)) write(u6,'(2a,i4,5x,a,f4.2,a)') 'GridName= ','VB orbital',iActOrb,' (',VBocc,')'
      end if
    end if
  end if
end do

if (isSphere) then
  !VV FIXME: no packing.
  write(line,'(a,i4)') 'Title= ',-1
  call PrintLine(LuVal,line,12,.true.)
  if (iBinary == 0) then
    do j=1,mCoor
      write(LuVal,'(E18.12)') SphrDist(j)
    end do
  else
    write(LuVal) (SphrDist(j),j=1,mCoor)
  end if
end if

if (isColor) then
  !VV FIXME: no packing.
  write(line,'(a,i4)') 'Title= ',-10
  call PrintLine(LuVal,line,12,.true.)
  if (iBinary == 0) then
    do j=1,mCoor
      write(LuVal,'(E18.12)') SphrColor(j)
    end do
  else
    write(LuVal) (SphrColor(j),j=1,mCoor)
  end if
end if

if (isDensity) then
  call outmo(0,2,WipMO,WipOcc,WipOut,mCoor,nMOs)
  do j=1,mCoor
    dNorm = dNorm+WipOut(j)
  end do
  !****
  !write(u6,*) ' mCoor=',mCoor
  !****
  if ((.not. isLine) .and. (.not. isLuscus)) then
    write(line,'(a,i4)') 'Title= ',0
    call PrintLine(LuVal,line,12,.true.)
  end if
  !if (isMOPack) then
  !  write(u6,*) 'pack code'
  !  call PackBlock(WipOut,iWipPBlock,mCoor,xLimits,iYDelta)
  !  write(line,9000) 0,(xLimits(j),j=1,4),(iYDelta(j),j=1,3)
  !  call PrintLine(LuVal,line,73,.false.)
  !
  !  if (iBinary /= 0) then
  !    call IArrToChar(iWipPBlock,cMoBlock,mCoor)
  !    !vv!!!!!!!
  !    if (isLuscus) then
  !      RC = C_WRITE(LID,CMOBLOCK,(mCoor*nBytesPackedVal)*RtoB) !!!!!!!!!!!!!!!!!!!!check mCoor*nBytesPackedVal
  !      if (RC == 0) THEN
  !        write(u6,*) 'error in writing luscus file!'
  !        call Abend()
  !      end if
  !    else
  !      write(LuVal) (cMoBlock(j),j=1,mCoor*nBytesPackedVal)
  !    end if
  !  else
  !    if (isLuscus) then
  !      RC = C_WRITE(LID,IWIPPBLOCK,(mCoor)*RtoB) !!!!!!!!!!!!!!!!!!!!check mCoor*nBytesPackedVal
  !      if (RC == 0) then
  !        write(u6,*) 'error in writing luscus file!'
  !        call Abend()
  !      end if
  !    else
  !      write(LuVal,'(I5)') (iWipPBlock(j), j=1,mCoor)
  !    end if
  !  endif
  !else
  if (isLuscus) then
    call dump_lusc(LID,WipOut,mCoor)
  end if

  if (iBinary == 1) then
    ! packing late
    if (isCutOff) then
      call mma_allocate(CMP,mCoor,label='TMP')
      iii = 0
      do ii=1,mCoor
        if (iWipCutOff(ii) == 1) then
          iii = iii+1
          CMP(iii) = WipOut(ii)
        end if
      end do

      if (isLuscus) then
        !!!!!!!!!!!!!!!!!!!!check iii-1
        RC = C_WRITE(LID,CMP,(III-1)*RtoB)
        if (RC == 0) then
          write(u6,*) 'error in writing luscus file!'
          call Abend()
        end if
      else
        write(LuVal) CMP(1:iii)
      end if
      call mma_deallocate(CMP)
    else
      ! no cut off
      !if (isLuscus) then
      !  call dump_lusc(LID,WipOut,mCoor)
      !  write(u6,*) 'here'
      !  RC = C_WRITE(LID,WIPOUT,MCOOR*RtoB) !!!!!!!!!!!!!!!!!!!!check MCOOR
      !  if (RC == 0) then
      !    write(u6,*) 'error in writing luscus file!'
      !    call Abend()
      !  end if
      !else
      write(LuVal) (WipOut(j),j=1,mCoor)
      !end if
    end if
  else

    if (isLine) then
      do j=1,mCoor
        WLine(1,j) = WipOut(j)
      end do
    else
      if (isDebug) then
        ! extra output -
        if (.not. isLuscus) write(LuVal,'(E18.12,3f8.4)') (WipOut(j),WCoor(1,j),WCoor(2,j),Wcoor(3,j),j=1,mCoor)
      else
        ! normal output - just numbers
        if (isCutOff) then
          do j=1,mCoor
            if (iWipCutOff(j) == 1) write(LuVal,'(E18.12)') WipOut(j)
          end do
        else

          if (.not. isLuscus) write(LuVal,'(E18.12)') (WipOut(j),j=1,mCoor)
        end if
      end if

    end if
    !GG This is only for testing CASDFT functional. It will be restore.
    !GG write(LuVal,'(E10.4)') (WipOut(j),j=1,mCoor)
  end if
end if
!end if

if (isLine .and. (.not. isLuscus)) then
  do i=1,mCoor
    write(LuVal,'(3F10.6,22E20.12)') (WCoor(j,i),j=1,3),(WLine(j,i),j=1,nLine)
  end do
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(isMOPack)
  call Unused_integer_array(iWipPBlock)
  call Unused_character_array(cMoBlock)
  call Unused_integer(nBytesPackedVal)
end if

end subroutine DumpM2Msi
