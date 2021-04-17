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
       subroutine zasun_zr  (i1,length,
     &                       valn,jn,kn,ln)
c
c     this routine write one block of 3-indices and appropriate
c     values of integrals into an opened TEMP-file
c
c     i1 - number of pivot index (I)
c     length - number of valid integrals in block (I)
c     this routine has also jn,kn,ln,valn
c     and stattemp and tmpnam as inputs, but they are
c     transpotred through commons  in reorg.fh
c
       implicit real*8 (a-h,o-z)
       integer length,i1
#include "reorg.fh"

#include "SysDef.fh"
       real*8 valn(1:nsize,1:mbas)
       integer jn(1:nsize,1:mbas)
       integer kn(1:nsize,1:mbas)
       integer ln(1:nsize,1:mbas)
c
c     help variable
c
       integer m2 ! ,iRec
c
       integer jkl(1:nsize)
       integer constj
       parameter (constj=1048576)
       integer constk
       parameter (constk=1024)
       integer f_iostat
       logical is_error
c
c*     pack indexes
c
       do m2=1,length
          jkl(m2)=ln(m2,i1)+constj*jn(m2,i1)
          jkl(m2)=jkl(m2)+constk*kn(m2,i1)
       end do
c
c*     open corresponding TEMP file in corresponding form
c
       if (iokey.eq.1) then
c
c      Fortran IO
c
       if (stattemp(i1).eq.0) then
c        file will be opened first time, it must be opened
c        whith the pointer at then first possition
         call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
c         open (unit=lunpublic,
c     &         file=tmpnam(i1),
c     &         form='unformatted',
c     &         status='unknown')
         stattemp(i1)=1
c
       else
c        file was alredy used in expansion of this block, it must
c        be opened whith the pointer at the end of the file
c@#ifdef _DECAXP_
       call molcas_open_ext2(lunpublic,tmpnam(i1),'append',
     &   'unformatted',f_iostat,.false.,1,'unknown',is_error)
cvv         open (unit=lunpublic,
cvv     &         file=tmpnam(i1),
cvv     &         form='unformatted',
cvv     &         status='unknown',
cvv     &         access='append')
c@#else
c@      call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
c         open (unit=lunpublic,
c     &         file=tmpnam(i1),
c     &         form='unformatted',
c     &         status='unknown')
c@       Do iRec = 1,nrectemp(i1)
c@         Read (lunpublic) m2
c@       End Do
c@#endif
c
       end if
c
       write (lunpublic) (valn(i,i1),i=1,length),
     &                   (jkl(i),i=1,length)
       close (lunpublic)
c
       else
c
c      MOLCAS IO
c
       call daname (lunpublic,tmpnam(i1))
       call ddafile (lunpublic,1,valn(1,i1),length,stattemp(i1))
       call idafile (lunpublic,1,jkl,       length,stattemp(i1))
       call daclos (lunpublic)
c
       end if
c
       nrectemp(i1)=nrectemp(i1)+1
       lrectemp(i1)=length
c
       return
       end
