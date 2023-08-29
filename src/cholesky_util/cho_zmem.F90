!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_ZMem(irc,l_Z,NVT,l_NVT,DoPrint,DoCheck)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: compute dimension of Z vector array and print it to
!          output (if requested through DoPrint). Check if there
!          is sufficient memory at this point (if DoCheck).
!
! Return codes:
!
! irc=-1 : Input error (debug only)
! irc=0  : All ok
! irc=1  : Negative length of Z vector array (could be integer
!          overflow)
! irc=999: Insufficient memory for Z vectors (only if DoCheck)

use Cholesky, only: LuPri, nSym
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc, l_Z
integer(kind=iwp), intent(in) :: l_NVT, NVT(l_NVT)
logical(kind=iwp), intent(in) :: DoPrint, DoCheck
integer(kind=iwp) :: iSym, l_Mx
real(kind=wp) :: Byte, Word(8), xl_Z
character(len=2) :: Unt
#if !defined (_I8_) || defined (_DEBUGPRINT_)
character(len=*), parameter :: SecNam = 'Cho_ZMem'
#endif

#ifdef _DEBUGPRINT_
if (l_NVT < nSym) then
  irc = -1
  l_Z = -999999
  return
end if
#endif

irc = 0

xl_Z = Zero
do iSym=1,nSym
  Word(iSym) = real(NVT(iSym),kind=wp)*(real(NVT(iSym),kind=wp)+One)*Half
  xl_Z = xl_Z+Word(iSym)
end do
l_Z = int(xl_Z)

if (DoPrint) then
  call Cho_Head('Z Vector Storage Requirements','-',80,LuPri)
  write(LuPri,*)
  do iSym=1,nSym
    call Cho_RWord2Byte(Word(iSym),Byte,Unt)
    write(LuPri,'(A,I2,A,I8,A,F8.3,1X,A,A)') 'Symmetry',iSym,':   ',int(Word(iSym)),' words (',Byte,Unt,')'
  end do
  write(LuPri,'(A)') '------------------------------------------'
  call Cho_RWord2Byte(xl_Z,Byte,Unt)
  write(LuPri,'(A,I8,A,F8.3,1X,A,A)') 'Total:        ',l_Z,' words (',Byte,Unt,')'
end if

#if !defined (_I8_) || defined (_DEBUGPRINT_)
if (l_Z < 0) then
  write(Lupri,'(A,A)') SecNam,': dimension of Z vector array is negative!'
  write(Lupri,'(A,I8)') 'l_Z=',l_Z
  if (xl_Z > Zero) then
    write(LuPri,'(A)') 'This seems to be an integer overflow!'
    call Cho_RWord2Byte(xl_Z,Byte,Unt)
    write(LuPri,'(A,1P,D15.6,A,D15.6,1X,A,A)') 'In double precision, xl_Z=',xl_Z,' words (',Byte,Unt,')'
  end if
  irc = 1
  return
end if
#endif

if (DoCheck) then
  call mma_maxDBLE(l_Mx)
  if (l_Z > l_Mx) then
    irc = 999
    return
  end if
end if

end subroutine Cho_ZMem
