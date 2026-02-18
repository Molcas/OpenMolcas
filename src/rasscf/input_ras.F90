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

module Input_RAS

use Definitions, only: iwp

implicit none
private

! Logical unit number for reading input:
integer(kind=iwp) :: LuInput
! Used for input processing
integer(kind=iwp), parameter :: NKeys = 149
! Logical flags, to check whether a keyword has been used in the input:
logical(kind=iwp) :: KeyFlags(0:NKeys)
! Actual keywords
character(len=4), parameter :: CMD(nKeys) = ['ALTE','ATOM','CHAR','CHOI','CHOL','CIMX','CION','CIRE','CIRO','CISE','CLEA','CORE', &
                                             'DELE','END ','FILE','FROZ','HOME','INAC','INPO','IPHN','ITER','JOBI','FUNC','LEVS', &
                                             'LINE','LOWD','LOWM','LUMO','MAXO','NACT','NEWI','NONE','NOQU','OPTO','ORBA','ORBL', &
                                             'ORBO','ORDE','OUTO','OUTP','PRIN','PROR','PRSP','PRWF','QUNE','RAS1','RAS2','RAS3', &
                                             'GASS','RFPE','CIRF','RFRO','RLXR','RASS','SDAV','SPIN','SUPS','SXDA','SYMM','THRS', &
                                             'TIGH','TITL','TYPE','VB  ','EXPE','SPLI','NUSP','ENSP','PESP','FOSP','MDRL','OFEM', &
                                             'FTHA','DFMD','BKAP','ALPH','FARO','DMRG','3RDM','NECI','TOTA','TIME','NMCY','CALC', &
                                             'RDMS','REAL','DEFI','DIAG','EMBD','SOCC','RGIN','PRSD','FCID','NOCA','SAVE','EXPA', &
                                             'H5OR','H5CI','HEXS','HEUR','DMPO','NEVP','HFOC','DAVT','CHRE','CHBL','MXSW','NOIS', &
                                             'DMRE','MXCA','DEXS','HROO','TDM ','DFCF','NKEE','REOR','TRIA','POPS','SEMI','MEMO', &
                                             'IVO ','CRPR','RDML','ORTH','CCCI','ROST','XMSI','CMSI','CMMA','CMMI','CMTH','GUGA', &
                                             'CMSS','CMSO','PERI','SSCR','MCM7','WRMA','DICE','STOC','EPSI','SAMP','DITE','DIRE', &
                                             'DIOC','PPT2','NDPT','RGRA','STAV']

public :: CMD, Key, KeyFlags, LuInput, NKeys, SetKey

contains

pure function Key(Keyword)

  logical(kind=iwp) :: Key
  character(len=*), intent(in) :: Keyword
  integer(kind=iwp) :: i

  do i=1,nKeys
    if (Keyword == Cmd(i)) then
      Key = KeyFlags(i)
      return
    end if
  end do
  Key = .false.

end function Key

subroutine SetKey(Keyword,Val)

  use Definitions, only: u6

  character(len=*), intent(in) :: Keyword
  logical(kind=iwp), intent(in) :: Val
  integer(kind=iwp) :: i

  do i=1,nKeys
    if (Keyword == Cmd(i)) then
      KeyFlags(i) = Val
      return
    end if
  end do

  write(u6,*) 'Keyword "',trim(Keyword),'" is unknown.'
  call Abend()

end subroutine SetKey

end module Input_RAS
