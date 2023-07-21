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

subroutine CHO_MCA_GETKEY(LUNIT,OPTION,LOPTION,NOPTION,IDKEY,LUPRI)
!
! Purpose: get next keyword and convert it to internal IDKEY.

implicit real*8(a-h,o-z)
character*(*) OPTION(NOPTION)  ! <-- character*(loption)
character*14 SECNAM
parameter(SECNAM='CHO_MCA_GETKEY')
parameter(LBLINE=180,LKWORD=180)
character*(LBLINE) BLINE
character*(LKWORD) KWORD, KEY
character*(LKWORD) GET_LN
external GET_LN
character*1 COMMENT
parameter(COMMENT='*')
integer CHO_TABIND
external CHO_TABIND
logical USE_OBS, USE_ALI
parameter(NTABLE=58,LKEY=4,NEOINP=1)
character*(LKEY) TABLE(NTABLE)
character*(LKEY) EOINP(NEOINP)
parameter(NOBSOL=1,NALIAS=12)
character*(LKEY) OBSOL(NOBSOL)
character*(LKEY) ALIAS(NALIAS,2)
data TABLE/'THRC','PRIN','BUFF','THRD','DMP1','DMP2','SPAN','MINQ','MXSH','SCRE','NOSC','QUAL','THRN','WARN','TOON','CHEC','CHKA', &
           'RSTD','RSTC','RSTM','MAXQ','CHOM','REDM','CHKS','CHKM','ABSO','NOAB','TRCN','IOVE','REOR','HALT','FRAC','QFRA','MXSU', &
           'ADDR','IFCS','ONES','TWOS','NAIV','VBUF','DIAC','TSTS','SSCR','NOSS','SSTH','SSNO','1-CE','NO2-','NOPR','PRES','PRET', &
           'PARA','SIMP','SIMR','FAKE','TRUE','BLOC','IDLE'/
data EOINP/'ENDC'/
data OBSOL/'XXXX'/
data ALIAS/'PREC','THSI','THSU','STOP','MEMQ','IOMO','DYNA','1CEN','NO2C','THRP','1CCD','1C-C','THRC','DMP1','DMP2','HALT','QFRA', &
           'ADDR','VBUF','1-CE','NO2-','PRET','1-CE','1-CE'/
! Set flags for using obsolete/alias keywords:
parameter(USE_OBS=.false.,USE_ALI=.true.)

! Check that we are in sync with caller.
! --------------------------------------

if (NOPTION /= NTABLE) then
  write(LUPRI,*) SECNAM,': NOPTION = ',NOPTION,' NTABLE = ',NTABLE
  IDKEY = -5
  GO TO 2000
end if

! Other internal checks.
! ----------------------

if (LKEY > LKWORD) then
  write(LUPRI,*) SECNAM,': LKEY = ',LKEY,' LKWORD = ',LKWORD
  IDKEY = -5
  GO TO 2000
end if

! Set blank line.
! ---------------

do I=1,LBLINE
  BLINE(I:I) = ' '
end do

! Read keyword (the check for comment/blank line should be
! obsolete).
! --------------------------------------------------------

1 KEY = GET_LN(LUNIT)
KWORD = KEY
call UPCASE(KWORD)
KWORD = adjustl(KWORD)
do while ((KWORD(1:1) == COMMENT) .or. (KWORD == BLINE))
  KEY = GET_LN(LUNIT)
  KWORD = KEY
  call UPCASE(KWORD)
  KWORD = adjustl(KWORD)
end do

LAST = ICLAST(KWORD,len(KWORD))
do I=LAST+1,LKEY
  KWORD(I:I) = ' '
end do

! Check for obsolete keyword.
! ---------------------------

if (USE_OBS) then
  IOBSOL = CHO_TABIND(OBSOL,LKEY,NOBSOL,' ',0,0,KWORD(1:LKEY))
  if ((IOBSOL > 0) .and. (IOBSOL <= NOBSOL)) then
    write(LUPRI,*) '*** NOTICE: Cholesky keyword "',KWORD(1:LKEY),'" is obsolete and will be disregarded.'
    GO TO 1
  end if
end if

! Check for alias.
! ----------------

IALIAS = 0
if (USE_ALI) then
  IALIAS = CHO_TABIND(ALIAS(1,1),LKEY,NALIAS,' ',0,0,KWORD(1:LKEY))
  if ((IALIAS > 0) .and. (IALIAS <= NALIAS)) then
    KWORD(1:LKEY) = ALIAS(IALIAS,2)
  else
    IALIAS = 0
  end if
end if

! Table lookup.
! -------------

IDKEY = CHO_TABIND(TABLE,LKEY,NTABLE,EOINP,LKEY,NEOINP,KWORD(1:LKEY))
if (IDKEY == -1) then
  write(LUPRI,*) SECNAM,': keyword not recognized:'
  if (IALIAS > 0) then
    write(LUPRI,*) 'Internal  key: ',KWORD(1:LAST),' (significant part: ',KWORD(1:LKEY),')'
    write(LUPRI,*) 'Aliasing used: ',ALIAS(IALIAS,1),' <-> ',ALIAS(IALIAS,2)
  else
    write(LUPRI,*) 'Internal  key: ',KWORD(1:LAST),' (significant part: ',KWORD(1:LKEY),')'
  end if
  write(LUPRI,*)
  if (LOPTION > 0) then
    write(LUPRI,*) 'Available keywords and short explanations:'
    do IOPTION=1,NOPTION
      if (TABLE(IOPTION) /= 'XXXX') write(LUPRI,*) TABLE(IOPTION),': ',OPTION(IOPTION)
    end do
  else
    write(LUPRI,*) 'Available keywords:'
    do IOPTION=1,NOPTION
      if (TABLE(IOPTION) /= 'XXXX') write(LUPRI,*) TABLE(IOPTION)
    end do
  end if
  write(LUPRI,*)
end if

! Normal exit point.
! ------------------

2000 continue
return

end subroutine CHO_MCA_GETKEY
