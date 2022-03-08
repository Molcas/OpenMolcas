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

subroutine DumpM2Msi(iRun,Luval,LID,nShowMOs,isDensity,nMOs,GRef,Occ,MO,DOut,mCoor,iGauss,nInc,isMOPack,PBlock,cMoBlock, &
                     nBytesPackedVal,dnorm,Crypt,isTheOne,isLine,iBinary,isEner,iType,NZ,E,WLine,nLine,WCoor,iPrintCount,isDebug, &
                     isCutOff,iCutOff,isSphere,SphrDist,isColor,SphrColor,isLuscus,NBYTES,NINLINE)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

#include "intent.fh"

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: iRun, LuVal, LID, nShowMOs, nMOs, GRef(*), mCoor, iGauss, nInc, PBlock(*), &
                                 nBytesPackedVal, iBinary, iType(*), NZ(*), nLine, iCutOff(*), NBYTES, NINLINE
logical(kind=iwp), intent(in) :: isDensity, isMOPack, isTheOne, isLine, isEner, isDebug, isCutOff, isSphere, isColor, isLuscus
integer(kind=iwp), intent(inout) :: iPrintCount
real(kind=wp), intent(in) :: Occ(*), MO(*), E(*), WCoor(3,mCoor), SphrDist(mCoor), SphrColor(mCoor)
real(kind=wp), intent(_OUT_) :: DOut(*)
character, intent(in) :: cMoBlock(*)
real(kind=wp), intent(inout) :: dNorm, WLine(nLine,mCoor)
character(len=7), intent(in) :: Crypt
integer(kind=iwp) :: i, ib, ii, iii, iMOs, j, RC !, iActOrb, iYDelta(3)
real(kind=wp) :: DumArr(2) !, xLimits(4)
character :: bb
character(len=128) :: Line
real(kind=wp), allocatable :: CMP(:)
!character(len=20) :: formt
! test
!character(len=3) :: cint
!character, parameter :: cx(64) = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ@#'
! test end
#include "macros.fh"
unused_var(isMOPack)
unused_var(PBlock(1))
unused_var(cMoBlock(1))
unused_var(nBytesPackedVal)

!write(u6,*) 'entering DumpM2Msi'
if (irun > 100) write(u6,*) iGauss,nbytes,nInc,ninline
!iActOrb = 0
iPrintCount = iPrintCount+1
do i=1,nShowMOs-merge(1,0,isDensity)-merge(1,0,isSphere)-merge(1,0,isColor)
  iMOs = GRef(i)

  !if ((Occ(iMOs) > Zero) .and. (Occ(iMOs) < Two)) then
  !  iActOrb = iActOrb+1
  !  call outmo(0,1,MO,VBmat(1+(iActOrb-1)*nMOs),DOut,mCoor,nMOs)
  !else
  call outmo(iMOs,1,MO,dumArr,DOut,mCoor,nMOs)
  !end if

  if ((.not. isLine) .and. (.not. isLuscus)) then
    write(line,'(a,i4)') 'Title= ',iMOs
    call PrintLine(LuVal,line,12,.true.)
  end if
  if (isTheOne) then
    if ((.not. isLine) .and. (.not. isLuscus)) then
      write(LuVal,'(f18.12)') (DOut(j),j=1,mCoor)
    else
      if (i+1 <= nLine) then
        do j=1,mCoor
          WLine(i+1,j) = DOut(j)
        end do
      end if
    end if
  else

    !if (isMOPack) then
    !  !call PackBlock(DOut,PBlock,mCoor,xLimits,iYDelta)
    !  write(line,9000) 0,(xLimits(j),j=1,4),(iYDelta(j),j=1,3)
    !  call PrintLine(LuVal,line,73,.false.)
    !  9000 format ('BHeader=',I2,1X,(4(E10.4,1X),3(I5,1X)))
    !  if (iBinary /= 0) then
    !    !call IArrToChar(PBlock,cMoBlock,mCoor)
    !    !vv ! NOT CODED YET
    !    write(LuVal) (cMoBlock(j),j=1,mCoor*nBytesPackedVal)
    !  else
    !    write(LuVal,'(I5)') (PBlock(j), j=0,mCoor-1)
    !  end if
    !else
    if (iBinary == 1) then
      ! packing late
      if (isCutOff) then
        call mma_allocate(CMP,mCoor,label='TMP')
        iii = 0
        do ii=1,mCoor
          if (iCutOff(ii) == 1) then
            iii = iii+1
            CMP(iii) = DOut(ii)
          end if
        end do

        write(LuVal) CMP(1:iii)
        call mma_deallocate(CMP)
      else
        ! no cut off
        write(LuVal) (DOut(j),j=1,mCoor)

      end if
    else !iBinary
      !write(cint,'(i3.3)') i
      if (isDebug) then
        ! extended output -
        write(LuVal,'(E10.4,3f8.4)') (DOut(j),WCoor(1,j),WCoor(2,j),Wcoor(3,j),j=1,mCoor)
      else !isDebug
        ! normal output - just numbers
        if (isCutOff) then
          do j=1,mCoor
            if (iCutOff(j) == 1) write(LuVal,'(E10.4)') DOut(j)
          end do
        else if (isLUSCUS) then
          ! NOPACKING
          !if (.not. isMOPack) then
          call dump_lusc(LID,DOut,mCoor)
          ! debug dump of data
          !
          !  do j=1,mCoor,NINLINE
          !    num = mCoor-j
          !    if (num > NINLINE) num = NINLINE
          !    write(formt,'(A,I2,A,I2,A,A)') '(',NINLINE,'E',NBYTES,'.4',')'
          !    write(line,formt) (DOut(j-1+ij),ij=1,num)
          !    call printline(LID,line,NINLINE*NBYTES,.false.)
          !  end do
          !else
          !  ! PACKING NOT implemented
          !  !open(unit=38,file='testgr_'//cint//'.txt')
          !  do j=1,min(mcoor,10)
          !    iexpnt = int(log10(abs(DOut(j))))
          !    write(u6,*) 'num=',DOut(j),'exp = ',iexpnt
          !    dnum = (One+DOut(j)/Ten**iexpnt)/Two
          !
          !    in1 = int(dnum*64.0_wp)
          !    in2 = int((dnum-in1*64.0_wp)*4096.0_wp)
          !    iexpnt = iexpnt+50
          !    write(u6,'(1x,3(1x,e18.8),2x,3(1x,i3))') DOut(j),dnum,dexpnt,in1,in2,iexpnt
          !    if (iexpnt < 1) then
          !      iexpnt = 1
          !    else if (iexpnt > 64) then
          !      iexpnt = 64
          !    end if
          !    write(u6,'(1x,3(1x,e18.8),2x,3(1x,i3))') DOut(j),dnum,dexpnt,in1,in2,iexpnt
          !    !                                       cx(in1:in1),cx(in2:in2),cx(iexpnt:iexpnt)
          !    !write(u6,*) '-----------------------'
          !  end do
          !  !close(38)
          !  !xxxmin = huge(xxxmin)
          !  !xxxmax = -huge(xxxmax)
          !  !do j=1,mCoor
          !  !  if (DOut(j) > xxxmax) xxxmax = DOut(j)
          !  !  if (DOut(j) < xxxmin) xxxmin = DOut(j)
          !  !end do
          !  !write(u6,*) 'test min/max',xxxmin,xxxmax
          !
          !end if !isMOPack
        else !isCutOff
          ! writing of data
          write(LuVal,'(E10.4)') (DOut(j),j=1,mCoor)
        end if !isCutOff, isLuscus
      end if !isDebug
    end if !iBinary
    !end if !isMOPack

    j = GRef(i)

    if (isEner) then
      !if ((Occ(j) > Zero) .and. (Occ(j) < Two)) then
      !  !iActOrb = iActOrb+1
      !  if ((iRun == 1) .and. (iPrintCount == 1)) write(u6,'(2a,i4,5x,a,f4.2,a)') 'GridName= ','VB orbital',iActOrb,' (',VBocc,')'
      !else
      ib = iType(j)
      bb = ' '
      if ((ib > 0) .and. (ib < 8)) bb = Crypt(ib:ib)
      if ((iRun == 1) .and. (iPrintCount == 1)) &
        write(u6,'(a,i2,i5,f12.4," (",f4.2,") ",a)') 'GridName= ',NZ(j),NZ(j+nMOs),E(j),Occ(j),bb
      !end if
    else
      !if ((Occ(j) > Zero) .and. (Occ(j) < Two)) then
      !  !iActOrb = iActOrb+1
      !  if ((iRun == 1) .and. (iPrintCount == 1)) write(u6,'(2a,i4,5x,a,f4.2,a)') 'GridName= ','VB orbital',iActOrb,' (',VBocc,')'
      !else
      ib = iType(j)
      bb = ' '
      if ((ib > 0) .and. (ib < 8)) bb = Crypt(ib:ib)
      if ((iRun == 1) .and. (iPrintCount == 1)) write(u6,'(a,i2,i5," (",f8.6,") ",a)') 'GridName= ',NZ(j),NZ(j+nMOs),Occ(j),bb
      !end if
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
  call outmo(0,2,MO,Occ,DOut,mCoor,nMOs)
  do j=1,mCoor
    dNorm = dNorm+DOut(j)
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
  !  call PackBlock(DOut,PBlock,mCoor,xLimits,iYDelta)
  !  write(line,9000) 0,(xLimits(j),j=1,4),(iYDelta(j),j=1,3)
  !  call PrintLine(LuVal,line,73,.false.)
  !
  !  if (iBinary /= 0) then
  !    call IArrToChar(PBlock,cMoBlock,mCoor)
  !    !vv!!!!!!!
  !    if (isLuscus) then
  !      RC = C_WRITE_WRAPPER(LID,CMOBLOCK,(mCoor*nBytesPackedVal)*RtoB) !!!!!!!!!!!!!!!!!!!!check mCoor*nBytesPackedVal
  !      if (RC == 0) THEN
  !        write(u6,*) 'error in writing luscus file!'
  !        call Abend()
  !      end if
  !    else
  !      write(LuVal) (cMoBlock(j),j=1,mCoor*nBytesPackedVal)
  !    end if
  !  else
  !    if (isLuscus) then
  !      RC = C_WRITE_WRAPPER(LID,PBLOCK,(mCoor)*RtoB) !!!!!!!!!!!!!!!!!!!!check mCoor*nBytesPackedVal
  !      if (RC == 0) then
  !        write(u6,*) 'error in writing luscus file!'
  !        call Abend()
  !      end if
  !    else
  !      write(LuVal,'(I5)') (PBlock(j), j=1,mCoor)
  !    end if
  !  endif
  !else
  if (isLuscus) then
    call dump_lusc(LID,DOut,mCoor)
  end if

  if (iBinary == 1) then
    ! packing late
    if (isCutOff) then
      call mma_allocate(CMP,mCoor,label='TMP')
      iii = 0
      do ii=1,mCoor
        if (iCutOff(ii) == 1) then
          iii = iii+1
          CMP(iii) = DOut(ii)
        end if
      end do

      if (isLuscus) then
        !!!!!!!!!!!!!!!!!!!!check iii-1
        RC = C_WRITE_WRAPPER(LID,CMP,(III-1)*RtoB)
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
      !  call dump_lusc(LID,DOut,mCoor)
      !  write(u6,*) 'here'
      !  RC = C_WRITE_WRAPPER(LID,DOUT,MCOOR*RtoB) !!!!!!!!!!!!!!!!!!!!check MCOOR
      !  if (RC == 0) then
      !    write(u6,*) 'error in writing luscus file!'
      !    call Abend()
      !  end if
      !else
      write(LuVal) (DOut(j),j=1,mCoor)
      !end if
    end if
  else

    if (isLine) then
      do j=1,mCoor
        WLine(1,j) = DOut(j)
      end do
    else
      if (isDebug) then
        ! extra output -
        if (.not. isLuscus) write(LuVal,'(E18.12,3f8.4)') (DOut(j),WCoor(1,j),WCoor(2,j),Wcoor(3,j),j=1,mCoor)
      else
        ! normal output - just numbers
        if (isCutOff) then
          do j=1,mCoor
            if (iCutOff(j) == 1) write(LuVal,'(E18.12)') DOut(j)
          end do
        else

          if (.not. isLuscus) write(LuVal,'(E18.12)') (DOut(j),j=1,mCoor)
        end if
      end if

    end if
    !GG This is only for testing CASDFT functional. It will be restore.
    !GG write(LuVal,'(E10.4)') (DOut(j),j=1,mCoor)
  end if
end if
!end if

if (isLine .and. (.not. isLuscus)) then
  do i=1,mCoor
    write(LuVal,'(3F10.6,22E20.12)') (WCoor(j,i),j=1,3),(WLine(j,i),j=1,nLine)
  end do
end if

return

contains

function c_write_wrapper(FileDescriptor,Buffer,nBytes)

  use, intrinsic :: iso_c_binding, only: c_loc

  integer(kind=iwp) :: c_write_wrapper
  integer(kind=iwp), intent(in) :: FileDescriptor, nBytes
  real(kind=wp), intent(in), target :: Buffer(*)
  interface
    function c_write(FileDescriptor,Buffer,nBytes) bind(C,name='c_write_')
      use, intrinsic :: iso_c_binding, only: c_ptr
      use Definitions, only: MOLCAS_C_INT
      integer(kind=MOLCAS_C_INT) :: c_write
      integer(kind=MOLCAS_C_INT) :: FileDescriptor, nBytes
      type(c_ptr), value :: Buffer
    end function c_write
  end interface

  c_write_wrapper = c_write(FileDescriptor,c_loc(Buffer(1)),nBytes)

end function c_write_wrapper

end subroutine DumpM2Msi
