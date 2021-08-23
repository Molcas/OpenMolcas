************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine filemanager (request,lun,rc)
c
c     request - specification of service (Input)
c     1 - open temporary file and give lun (minfiles-maxfiles)
c     2 - rewind file lun
c     3 - close file lun (and delete it if it is temporrary one)
c     4 - open file with given lun (used for fixed files 10-minfiles)
c     5 - close without deleting (for any file)
c     lun     - Logical unit number (Input, Output)
c     rc      - return (error) code (Output)
c
c     This routine is a manager of temporarry disk files, used durring
c     the calculations.
c     if request=1 it finds free lun number if it is possible
c     and open file as unformatted
c     if request=2 it rewinds lun file if it is opened
c     if request=3 if close given lun file if it is opened
c
#include "filemgr.fh"
#include "ccsd1.fh"
c
       integer lun,request,rc
c
c     help variables
c
       integer nhelp,mhelp,ierr,idum(1)
       logical is_error
c
       rc=0
c
c
       if (request.eq.1) then
c
cI    open new file
c
cI.1  look for lowest free lun
c
       do 100 nhelp=minfiles,maxfiles
        mhelp = nhelp
        if (filestatus(nhelp).eq.0) goto 101
 100    continue
c
cI.2  RC=1 : there is not enough Temporarry files allowed
       rc=1
       return
c
cI.3  open file
c
 101    lun=mhelp
c
       if (iokey.eq.1) then
c      Fortran IO
       call molcas_open_ext2(lun,filename(lun),'sequential',
     &  'unformatted',
     &      ierr,.false.,1,'unknown',is_error)
c       open (unit=lun,
c     &       file=filename(lun),
c     &       status='unknown',
c     &       form='unformatted',
c     &       iostat=ierr)
c
       else
c      MOLCAS IO
       call daname (lun,filename(lun))
       daddr(lun)=0
       end if
c
       filestatus(lun)=1
c
c
       else if (request.eq.2) then
c
cII   rewind given file
c
cII.1 test lun, if it is not too large
c
       if ((lun.lt.10).or.(lun.gt.maxfiles)) then
c     RC=2 : lun is out of range
       rc=2
       return
       end if
c
cII.2 test, if lun is opened
c
       if (filestatus(lun).ne.1) then
c     RC=3 : file with logical unit lun is not opened, it can't be rewined
       rc=3
       return
       end if
c
cII.3 rewind lun file
c
       if (iokey.eq.1) then
c      Fortran IO
       rewind (lun)
c
       else
c      MOLCAS IO
       nhelp=0
       idum(1)=nhelp
       call idafile (lun,5,idum,1,daddr(lun))
       end if
c
c
       else if (request.eq.3) then
c
cIII  close given file + scratch if it is Temp
c
cIII.1test lun, if it is not too large
c
       if ((lun.lt.10).or.(lun.gt.maxfiles)) then
c     RC=4 : lun is out of range
       rc=4
       return
       end if
c
cIII.2test, if lun is opened
c
       if (filestatus(lun).ne.1) then
c     RC=5 : file with logical unit lun is not opened,it can't be closed
       rc=5
       return
       end if
c
cIII.3 close ; scratch file, if it was a temporarry one
c
       if (iokey.eq.1) then
c      Fortran IO
         if (lun.ge.minfiles) then
c        close and scratch
         close (lun)
         call molcas_open(lun,filename(lun))
c         open (unit=lun,file=filename(lun))
         write (lun,*) ' File scratched'
         close (lun)
c        call sqname (lun,filename'lun')
c        call sqeras (lun)
         else
c        close only
         close (lun)
         end if
c
       else
c      MOLCAS IO
         if (lun.ge.minfiles) then
c        close and scratch
         call daeras (lun)
         else
c        close only
         call daclos (lun)
         end if
c
       end if
c
       filestatus(lun)=0
c
c
       else if (request.eq.4) then
c
cIV   open file with given lun
c
c
cIV.1 test lun, if it is not too large
c
       if ((lun.lt.10).or.(lun.gt.maxfiles)) then
c     RC=6 : lun is out of range
       rc=6
       return
       end if
c
cIV.2  test if file is not already opened
       if (filestatus(lun).eq.1) then
c     RC=7 : file is already opened
       rc=7
       return
       end if
c
cIV.3 open file
c
       if (iokey.eq.1) then
c      Fortran IO
       call molcas_open_ext2(lun,filename(lun),
     &     'sequential','unformatted',
     &      ierr,.false.,1,'unknown',is_error)
c       open (unit=lun,
c     &       file=filename(lun),
c     &       status='unknown',
c     &       form='unformatted',
c     &       iostat=ierr)
c
       else
c      MOLCAS IO
       call daname (lun,filename(lun))
       daddr(lun)=0
       end if
c
       filestatus(lun)=1
c
       else if (request.eq.5) then
c
cV     close file with given lun (without deleting)
c
cV.1   test lun, if it is not too large
c
       if ((lun.lt.10).or.(lun.gt.maxfiles)) then
c     RC=8 : lun is out of range
       rc=8
       return
       end if
c
cV.2  test, if lun is opened
c
       if (filestatus(lun).ne.1) then
c     RC=9 : file with logical unit lun is not opened,it can't be closed
       rc=9
       return
       end if
c
cV.3  close
c
       if (iokey.eq.1) then
c      Fortran IO
       close (lun)
c
       else
c      MOLCAS IO
       call daclos (lun)
       end if
c
       filestatus(lun)=0
c
c
       else
c     RC=10 : incorect value of request
       rc=10
       return
       end if
c
       return
       end
