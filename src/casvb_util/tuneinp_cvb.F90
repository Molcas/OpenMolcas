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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine tuneinp_cvb()

use casvb_global, only: alftol, cnrmtol, delopth1, delopth2, dfx, dfxmin, dfxtol, dx, eigwrngtol, endwhenclose, exp12tol, follow, &
                        grd, grdwrngtol, hhaccfac, hhrejfac, hhstart, hhtol, mxdav, nopth1, nopth2, nortiter, orththr, resthr, &
                        safety, scalesmall, sgn, signtol, singul, zzacclim, zzmax, zzmin, zzrejmax, zzrejmin

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iaux(1), istr, istr2, j, nread
real(kind=wp) :: daux(1)
integer(kind=iwp), parameter :: ncmp = 8, nstrin = 37
character(len=*), parameter :: string(nstrin) = ['CNRMTOL ','SAFETY  ','SIGNTOL ','ALFTOL  ','DFXTOL  ','EXP12TOL','GRDWRNGT', &
                                                 'EIGWRNG ','SINGUL  ','DFX     ','SIGN    ','ZZMAX   ','ZZMIN   ','DX      ', &
                                                 'GRD     ','NOPTH1  ','NOPTH2  ','DELOPTH1','DELOPTH2','HHREJFAC','HHACCFAC', &
                                                 'ZZACCLIM','HHTOL   ','DFXMIN  ','ZZREJMIN','ZZREJMAX','SCALESMA','HHSTART ', &
                                                 'RESTHR  ','NORTITER','ORTHTHR ','FOLLOW  ','MXDAV   ','LASTUPD ','ENDIFCLO', &
                                                 'ENDTUNE ','END     '], &
                               true(1) = ['T']

do
  call fstring_cvb(string,nstrin,istr,ncmp,2)
  if ((istr == 36) .or. (istr == 37)) then
    !'ENDTUNE' or 'END'
    istr = 0
  end if
  if (istr == 1) then
    ! 'CNRMTOL'
    call real_cvb(daux,1,nread,1)
    cnrmtol = daux(1)
  else if (istr == 2) then
    ! 'SAFETY'
    call real_cvb(daux,1,nread,1)
    safety = daux(1)
  else if (istr == 3) then
    ! 'SIGNTOL'
    call real_cvb(daux,1,nread,1)
    signtol = daux(1)
  else if (istr == 4) then
    ! 'ALFTOL'
    call real_cvb(daux,1,nread,1)
    alftol = daux(1)
  else if (istr == 5) then
    ! 'DFXTOL'
    call real_cvb(daux,1,nread,1)
    dfxtol = daux(1)
  else if (istr == 6) then
    ! 'EXP12TOL'
    call real_cvb(daux,1,nread,1)
    exp12tol = daux(1)
  else if (istr == 7) then
    ! 'GRDWRNGT'
    call real_cvb(daux,1,nread,1)
    grdwrngtol = daux(1)
  else if (istr == 8) then
    ! 'EIGWRNG '
    call real_cvb(daux,1,nread,1)
    eigwrngtol = daux(1)
  else if (istr == 9) then
    ! 'SINGUL'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 3)) then
      write(u6,*) ' Illegal I index in SINGUL :',i
      call abend_cvb()
    end if
    call real_cvb(singul(i),1,nread,1)
  else if (istr == 10) then
    ! 'DFX'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 6)) then
      write(u6,*) ' Illegal I index in DFX :',i
      call abend_cvb()
    end if
    call real_cvb(dfx(i),1,nread,1)
  else if (istr == 11) then
    ! 'SIGN'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 6)) then
      write(u6,*) ' Illegal I index in SIGN :',i
      call abend_cvb()
    end if
    call real_cvb(sgn(i),1,nread,1)
  else if (istr == 12) then
    ! 'ZZMAX'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 6)) then
      write(u6,*) ' Illegal I index in ZZMAX :',i
      call abend_cvb()
    end if
    call real_cvb(zzmax(i),1,nread,1)
  else if (istr == 13) then
    ! 'ZZMIN'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 6)) then
      write(u6,*) ' Illegal I index in ZZMIN :',i
      call abend_cvb()
    end if
    call real_cvb(zzmin(i),1,nread,1)
  else if (istr == 14) then
    ! 'DX'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 3)) then
      write(u6,*) ' Illegal I index in DX :',i
      call abend_cvb()
    end if
    call int_cvb(iaux,1,nread,1)
    j = iaux(1)
    if ((j < 1) .or. (j > 6)) then
      write(u6,*) ' Illegal J index in DX :',j
      call abend_cvb()
    end if
    call real_cvb(dx(i,j),1,nread,1)
  else if (istr == 15) then
    ! 'GRD'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 3)) then
      write(u6,*) ' Illegal I index in GRD :',i
      call abend_cvb()
    end if
    call int_cvb(iaux,1,nread,1)
    j = iaux(1)
    if ((j < 1) .or. (j > 6)) then
      write(u6,*) ' Illegal J index in GRD :',j
      call abend_cvb()
    end if
    call real_cvb(grd(i,j),1,nread,1)
  else if (istr == 16) then
    ! 'NOPTH1'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in NOPTH1 :',i
      call abend_cvb()
    end if
    call int_cvb(nopth1(i),1,nread,1)
  else if (istr == 17) then
    ! 'NOPTH2'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in NOPTH2 :',i
      call abend_cvb()
    end if
    call int_cvb(nopth2(i),1,nread,1)
  else if (istr == 18) then
    ! 'DELOPTH1'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in DELOPTH1 :',i
      call abend_cvb()
    end if
    call real_cvb(delopth1(i),1,nread,1)
  else if (istr == 19) then
    ! 'DELOPTH2'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in DELOPTH2 :',i
      call abend_cvb()
    end if
    call real_cvb(delopth2(i),1,nread,1)
  else if (istr == 20) then
    ! 'HHREJFAC'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in HHREJFAC :',i
      call abend_cvb()
    end if
    call real_cvb(hhrejfac(i),1,nread,1)
  else if (istr == 21) then
    ! 'HHACCFAC'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 5)) then
      write(u6,*) ' Illegal I index in HHACCFAC :',i
      call abend_cvb()
    end if
    call int_cvb(iaux,1,nread,1)
    j = iaux(1)
    if ((j < 1) .or. (j > 2)) then
      write(u6,*) ' Illegal J index in HHACCFAC :',j
      call abend_cvb()
    end if
    call real_cvb(hhaccfac(i,j),1,nread,1)
  else if (istr == 22) then
    ! 'ZZACCLIM'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 4)) then
      write(u6,*) ' Illegal I index in ZZACCLIM :',i
      call abend_cvb()
    end if
    call int_cvb(iaux,1,nread,1)
    j = iaux(1)
    if ((j < 1) .or. (j > 2)) then
      write(u6,*) ' Illegal J index in ZZACCLIM :',j
      call abend_cvb()
    end if
    call real_cvb(zzacclim(i,j),1,nread,1)
  else if (istr == 23) then
    ! 'HHTOL'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in HHTOL :',i
      call abend_cvb()
    end if
    call real_cvb(hhtol(i),1,nread,1)
  else if (istr == 24) then
    ! 'DFXMIN'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in DFXMIN :',i
      call abend_cvb()
    end if
    call real_cvb(dfxmin(i),1,nread,1)
  else if (istr == 25) then
    ! 'ZZREJMIN'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in ZZREJMIN :',i
      call abend_cvb()
    end if
    call real_cvb(zzrejmin(i),1,nread,1)
  else if (istr == 26) then
    ! 'ZZREJMAX'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in ZZREJMAX :',i
      call abend_cvb()
    end if
    call real_cvb(zzrejmax(i),1,nread,1)
  else if (istr == 27) then
    ! 'SCALESMALL'
    call int_cvb(iaux,1,nread,1)
    i = iaux(1)
    if ((i < 1) .or. (i > 2)) then
      write(u6,*) ' Illegal I index in SCALESMALL :',i
      call abend_cvb()
    end if
    call fstring_cvb(true,1,istr2,1,1)
    scalesmall(i) = (istr2 == 1)
  else if (istr == 28) then
    ! 'HHSTART'
    call real_cvb(daux,1,nread,1)
    hhstart = daux(1)
  else if (istr == 29) then
    ! 'RESTHR'
    call real_cvb(daux,1,nread,1)
    resthr = daux(1)
  else if (istr == 30) then
    ! 'NORTITER'
    call int_cvb(iaux,1,nread,1)
    nortiter = iaux(1)
  else if (istr == 31) then
    ! 'ORTHTHR'
    call real_cvb(daux,1,nread,1)
    orththr = daux(1)
  else if (istr == 32) then
    ! 'FOLLOW'
    call fstring_cvb(true,1,istr2,1,1)
    follow = (istr2 == 1)
  else if (istr == 33) then
    ! 'MXDAV'
    call int_cvb(iaux,1,nread,1)
    mxdav = iaux(1)
  !else if (istr == 34) then
  !  ! 'LASTUPD'
  !  call fstring_cvb(true,1,istr2,1,1)
  !  lastupd = (istr2 == 1)
  else if (istr == 35) then
    ! 'ENDIFCLOSE'
    call fstring_cvb(true,1,istr2,1,1)
    endwhenclose = (istr2 == 1)
  end if
  ! 'END', 'ENDTUNE' or unrecognized keyword -- end of input:
  if (istr == 0) exit
end do

return

end subroutine tuneinp_cvb
