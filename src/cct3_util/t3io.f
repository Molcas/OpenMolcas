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
c      cct3_getmediate
c      wrtmediate
c      cct3_getmap
c      cct3_wrtmap
c      rea
c      wri
c
c      all routines are analogous to those in ccsd
c
c      -------------------------------------------------------------
c
       subroutine cct3_getmediate (wrk,wrksize,
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
      call cct3_getmap (lun,poss0,lenght,mapd,mapi,rc1)
c
c2    read mediate in one block
c
       if (lenght.eq.0) then
c     RC=1 : there is nothing to read, lenght of mediate is 0
       rc=1
       return
       end if
c
       call cct3_rea (lun,lenght,wrk(poss0))
c
       return
       end
c
c     ----------------------------
c
       subroutine cct3_wrtmediate (wrk,wrksize,
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
      call cct3_wrtmap (lun,mapd,mapi,rc1)
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
       call cct3_wri (lun,lenght,wrk(poss0))
c
       return
       end
c
c     ----------------------------
c
       subroutine cct3_getmap (lun,poss0,lenght,mapd,mapi,rc)
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
#include "t31.fh"
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
       read (lun) mapd,mapi
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
       subroutine cct3_wrtmap (lun,mapd,mapi,rc)
c
c     this routine write required mapd and mapi to opened unformatted file
c     with number lun
c
c     lun   - Logical unit number of file, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Input)
c     mapi  - inverse map matrix corresponding to given mediate (Input)
c     rc    - return (error) code (Output)
c
#include "t31.fh"
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
       subroutine cct3_rea (lun,lenght,vector)
c
c     this routine read lenght-R8 numbers from opened unformatted file
c     with number lun form the given possition as one record
c
c     lun    - Logical unit number of file, where mediate is stored (Input)
c     lenght - # of R8 numbers to be read  (Input)
c     vector - space, where numbers are stored after reading  (Output)

c
#include "t31.fh"
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
       subroutine cct3_wri (lun,lenght,vector)
c
c     this routine write lenght-R8 numbers to opened unformatted file
c     with number lun at the given possition as one record
c
c     lun    - Logical unit number of file, where mediate will be stored (Input)
c     lenght - # of R8 numbers to be written  (Input)
c     vector - space, where numbers are stored  (Input)

c
#include "t31.fh"
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
