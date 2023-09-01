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

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LUNIT, LOPTION, NOPTION, LUPRI
character(len=LOPTION), intent(in) :: OPTION(NOPTION)
integer(kind=iwp), intent(out) :: IDKEY
integer(kind=iwp), parameter :: LKEY = 4, LKWORD = 180, NALIAS = 12, NEOINP = 1, NOBSOL = 1, NTABLE = 58
integer(kind=iwp) :: I, IALIAS, IOBSOL, IOPTION, LAST
character(len=LKWORD) :: KWORD, KEY
! Set flags for using obsolete/alias keywords:
logical(kind=iwp), parameter :: USE_ALI = .true., USE_OBS = .false.
character(len=*), parameter :: COMMENT = '*', SECNAM = 'CHO_MCA_GETKEY'
character(len=LKEY), parameter :: ALIAS(NALIAS,2) = reshape(['PREC','THSI','THSU','STOP','MEMQ','IOMO','DYNA','1CEN','NO2C', &
                                                             'THRP','1CCD','1C-C', &
                                                             'THRC','DMP1','DMP2','HALT','QFRA','ADDR','VBUF','1-CE','NO2-', &
                                                             'PRET','1-CE','1-CE'],[NALIAS,2]), &
                                  EOINP(NEOINP) = ['ENDC'], OBSOL(NOBSOL) = ['XXXX'], &
                                  TABLE(NTABLE) = ['THRC','PRIN','BUFF','THRD','DMP1','DMP2','SPAN','MINQ','MXSH','SCRE','NOSC', &
                                                   'QUAL','THRN','WARN','TOON','CHEC','CHKA','RSTD','RSTC','RSTM','MAXQ','CHOM', &
                                                   'REDM','CHKS','CHKM','ABSO','NOAB','TRCN','IOVE','REOR','HALT','FRAC','QFRA', &
                                                   'MXSU','ADDR','IFCS','ONES','TWOS','NAIV','VBUF','DIAC','TSTS','SSCR','NOSS', &
                                                   'SSTH','SSNO','1-CE','NO2-','NOPR','PRES','PRET','PARA','SIMP','SIMR','FAKE', &
                                                   'TRUE','BLOC','IDLE']
integer(kind=iwp), external :: ICLAST, CHO_TABIND
character(len=LKWORD), external :: GET_LN

! Check that we are in sync with caller.
! --------------------------------------

if (NOPTION /= NTABLE) then
  write(LUPRI,*) SECNAM,': NOPTION = ',NOPTION,' NTABLE = ',NTABLE
  IDKEY = -5
  return
end if

! Other internal checks.
! ----------------------

if (LKEY > LKWORD) then
  write(LUPRI,*) SECNAM,': LKEY = ',LKEY,' LKWORD = ',LKWORD
  IDKEY = -5
  return
end if

! Read keyword (the check for comment/blank line should be
! obsolete).
! --------------------------------------------------------

do
  KEY = GET_LN(LUNIT)
  KWORD = KEY
  call UPCASE(KWORD)
  KWORD = adjustl(KWORD)
  do while ((KWORD(1:1) == COMMENT) .or. (KWORD == ''))
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

  if (.not. USE_OBS) exit
  IOBSOL = CHO_TABIND(OBSOL,LKEY,NOBSOL,' ',0,0,KWORD(1:LKEY))
  if ((IOBSOL <= 0) .or. (IOBSOL > NOBSOL)) exit
  write(LUPRI,*) '*** NOTICE: Cholesky keyword "',KWORD(1:LKEY),'" is obsolete and will be disregarded.'
end do

! Check for alias.
! ----------------

IALIAS = 0
if (USE_ALI) then
  IALIAS = CHO_TABIND(ALIAS(:,1),LKEY,NALIAS,' ',0,0,KWORD(1:LKEY))
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

return

end subroutine CHO_MCA_GETKEY
