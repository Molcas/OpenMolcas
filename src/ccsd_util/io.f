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
c
c      this file contains following routines:
c
c      getmediate
c      wrtmediate
c      getmap
c      wrtmap
c      rea
c      wri
c      filemanager
c      mkfilemgrcom
c
c      -------------------------------------------------------------
c
       subroutine getmediate (wrk,wrksize,
     & lun,poss0,mapd,mapi,rc)
c
c     this routine read required mediate from opened unformatted file
c     with number lun, and place it statring with the poss0
c     it also reads mapd and mapi of the given mediade, and reconstruct
c     mapd to actual possitions
c
c     lun   - Logical unit number of file, where mediate is stored (Input)
c     poss0 - initial possition in WRK, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Output)
c     mapi  - inverse map matrix corresponding to given mediate (Output)
c     rc    - return (error) code (Output)
c
c     N.B.
c     all mediates are storred as follows
c     1 - mapd, mapi
c     2 - one record with complete mediate
c
#include "wrk.fh"
       integer lun,poss0,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer lenght,rc1
c
       rc=0
c
c1    read mapd
c
      call getmap (lun,poss0,lenght,mapd,mapi,rc1)
c
c2    read mediate in one block
c
       if (lenght.eq.0) then
c     RC=1 : there is nothing to read, lenght of mediate is 0
       rc=1
       return
       end if
c
       call rea (lun,lenght,wrk(poss0))
c
       return
       end
c
c     ----------------------------
c
       subroutine wrtmediate (wrk,wrksize,
     & lun,mapd,mapi,rc)
c
c     this routine write required mediate to opened unformatted file
c     with number lun
c     it also store mapd and mapi of the given mediade
c
c     lun   - Logical unit number of file, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Input)
c     mapi  - inverse map matrix corresponding to given mediate (Input)
c     rc    - return (error) code (Output)
c
c     N.B.
c     all mediates are storred as follows
c     1 - mapd, mapi
c     2 - one record with complete mediate
c
#include "wrk.fh"
       integer lun,rc,rc1
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer im,lenght,poss0
c
       rc=0
c
c1    write mapd
c
      call wrtmap (lun,mapd,mapi,rc1)
c
c2    calculate overall lenght
c
       lenght=0
c
       do 100 im=1,mapd(0,5)
       lenght=lenght+mapd(im,2)
 100    continue
c
c     write mediate in one block
c
       if (lenght.eq.0) then
c     RC=1 : there is nothing to write, lenght of mediate is 0
       rc=1
       return
       end if
c
       poss0=mapd(1,1)
       call wri (lun,lenght,wrk(poss0))
c
       return
       end
c
c     ----------------------------
c
       subroutine getmap (lun,poss0,lenght,mapd,mapi,rc)
c
c     this routine reads mapd and mapi of the given mediade
c     from lun and reconstruct mapd to actual possitions poss0
c
c     lun   - Logical unit number of file, where mediate is stored (Input)
c     poss0 - initial possition in WRK, where mediate will be stored (Input)
c     lenght- overall lenght of mediate (Output)
c     mapd  - direct map matrix corresponding to given mediate (Output)
c     mapi  - inverse map matrix corresponding to given mediate (Output)
c     rc    - return (error) code (Output)
c
c     N.B.
c     all mediates are storred as follows
c     1 - mapd, mapi
c     2 - one record with complete mediate
c
#include "filemgr.fh"
#include "ccsd1.fh"

#include "SysDef.fh"
c
       integer lun,poss0,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer poss,im,lenght
c
       rc=0
c
c1    read mapd
c
       if (iokey.eq.1) then
c      Fortran IO
       read (lun) ((mapd(i,j),i=0,512),j=1,6),
     & (((mapi(l,m,n),l=1,8),m=1,8),n=1,8)
c workaround for a bug in some Intel versions
#ifdef __INTEL_COMPILER
       else if (iokey.lt.0) then
       write(6,*) "this should never happen"
#endif
c
       else
c      MOLCAS IO
       call idafile (lun,2,mapd,513*6,daddr(lun))
       call idafile (lun,2,mapi,8*8*8,daddr(lun))
       end if
c
c2    change possitions in mapd to proper one and calculate overall lenght
c
       poss=poss0
       lenght=0
c
       do 100 im=1,mapd(0,5)
c
       mapd(im,1)=poss
       poss=poss+mapd(im,2)
       lenght=lenght+mapd(im,2)
c
 100    continue
c
       return
       end
c
c     ----------------------------
c
       subroutine wrtmap (lun,mapd,mapi,rc)
c
c     this routine write required mapd and mapi to opened unformatted file
c     with number lun
c
c     lun   - Logical unit number of file, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Input)
c     mapi  - inverse map matrix corresponding to given mediate (Input)
c     rc    - return (error) code (Output)
c
#include "filemgr.fh"
#include "ccsd1.fh"

#include "SysDef.fh"
c
       integer lun,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
       rc=0
c
c1    write mapd
c
       if (iokey.eq.1) then
c      Fortran IO
       write (lun) mapd,mapi
c
       else
c      MOLCAS IO
       call idafile (lun,1,mapd,513*6,daddr(lun))
       call idafile (lun,1,mapi,8*8*8,daddr(lun))
       end if
c
       return
       end
c
c     ----------------------------
c
       subroutine rea (lun,lenght,vector)
c
c     this routine read lenght-R8 numbers from opened unformatted file
c     with number lun form the given possition as one record
c
c     lun    - Logical unit number of file, where mediate is stored (Input)
c     lenght - # of R8 numbers to be read  (Input)
c     vector - space, where numbers are stored after reading  (Output)

c
#include "filemgr.fh"
#include "ccsd1.fh"

#include "SysDef.fh"
c
       integer lun,lenght
       real*8 vector(1:lenght)
c
       if (iokey.eq.1) then
c      Fortran IO
       read (lun) (vector(i),i=1,lenght)
c
       else
c      MOLCAS IO
       call ddafile (lun,2,vector,lenght,daddr(lun))
       end if
c
       return
       end
c
c     ----------------------------
c
       subroutine wri (lun,lenght,vector)
c
c     this routine write lenght-R8 numbers to opened unformatted file
c     with number lun at the given possition as one record
c
c     lun    - Logical unit number of file, where mediate will be stored (Input)
c     lenght - # of R8 numbers to be written  (Input)
c     vector - space, where numbers are stored  (Input)

c
#include "filemgr.fh"
#include "ccsd1.fh"

#include "SysDef.fh"
c
       integer lun,lenght
       real*8 vector(1:lenght)
c
       if (iokey.eq.1) then
c      Fortran IO
       write (lun) vector
c
       else
c      MOLCAS IO
       call ddafile (lun,1,vector,lenght,daddr(lun))
       end if
c
       return
       end
c
c     ----------------------------
c
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
       integer nhelp,mhelp,ierr
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
       call idafile (lun,5,nhelp,1,daddr(lun))
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
c
c     ----------------------------
c
       subroutine mkfilemgrcom
c
c     this routine make names of temp files and define initial values
c     of filestatus
c
#include "filemgr.fh"
#include "ccsd1.fh"
c
c     help variable
c
       integer nhelp
c
c1    def filestatus
c
       do 100 nhelp=10,maxfiles
       filestatus(nhelp)=0
 100    continue
c
c2    def disk addresses
c
       do 200 nhelp=10,maxfiles
       daddr(nhelp)=0
 200    continue
c
c3    def filenames
c
       filename(10)='INTAB'
       filename(11)='INTA1'
       filename(12)='INTA2'
       filename(13)='INTA3'
       filename(14)='INTA4'
       filename(15)='INTSTA'
       filename(16)=filerst
       filename(17)='Temp17'
       filename(18)='Temp18'
       filename(19)='Temp19'
       filename(20)='Temp20'
       filename(21)='Temp21'
       filename(22)='Temp22'
       filename(23)='Temp23'
       filename(24)='Temp24'
       filename(25)='Temp25'
       filename(26)='Temp26'
       filename(27)='Temp27'
       filename(28)='Temp28'
       filename(29)='Temp29'
       filename(30)='Temp30'
       filename(31)='Temp31'
       filename(32)='Temp32'
       filename(33)='Temp33'
       filename(34)='Temp34'
       filename(35)='Temp35'
       filename(36)='Temp36'
       filename(37)='Temp37'
       filename(38)='Temp38'
       filename(39)='Temp39'
       filename(40)='Temp40'
       filename(41)='Temp41'
       filename(42)='Temp42'
       filename(43)='Temp43'
       filename(44)='Temp44'
       filename(45)='Temp45'
       filename(46)='Temp46'
       filename(47)='Temp47'
       filename(48)='Temp48'
       filename(49)='Temp49'
       filename(50)='Temp50'
c

       return
       end
c
