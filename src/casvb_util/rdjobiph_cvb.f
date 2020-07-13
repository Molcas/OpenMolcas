************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine rdjobiph_cvb(fnjob)
      implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "rasdim.fh"
#include "jobiph_j.fh"

c Input parameters :
c ------------------
      character*(*) fnjob               ! subroutine arguments
c      integer mxtit,mxorb,mxroot,mxsym  ! rasdim.fh
c      integer itob,rtoi                 ! SysDef.fh
c Output parameters (in jobiph_j) :
c -------------------
c      integer iadr15_j(15)
c      integer nactel_j,ispin_j,nsym_j,lsym_j,nfro_j(mxsym),
c     >  nish_j(mxsym),nash_j(mxsym),ndel_j(mxsym),
c     >  nbas_j(mxsym)
c       character*(lenin8) name_j(mxorb)
c      integer nconf_j
c       character*2 header_j(72)
c       character*72 title_j(18)
c      real*8 potnuc_j
c      integer lroots_j,
c     >  nroots_j,iroot_j(mxroot),
c     >  nrs1_j(mxsym),nrs2_j(mxsym),nrs3_j(mxsym),
c     >  nhole1_j,nelec3_j,ipt2_j
c      real*8 weight_j(mxroot)
cw
cw    NA                   NORB                NNA              NTASH
cw    NTASQR             NTATRI              NTBAS             NTBSQR
cw    NTBTRI              NTISH             NTISQR             NTITRI
c Local/work parameters :
c -----------------------
c  local / work
      integer lujob
      integer idisk
      logical debug
      data debug/.false./
c      integer iwork        ! WrkSpc.fh

c Read the table of disk adresses :
      lujob=15
      call daname_cvb(lujob,fnjob)
      idisk=0
      call idafile(lujob,2,iadr15_j,15,idisk)
c Read the the system description :
      idisk=iadr15_j(1)
      call WR_RASSCF_Info(lujob,2,idisk,
     >                    nactel_j,ispin_j,nsym_j,lsym_j,nfro_j,
     >                    nish_j,nash_j,ndel_j,
     >                    nbas_j,mxsym,name_j,lenin8*mxorb,
     >                    nconf_j,header_j,144,
     >                    title_j,4*18*mxtit,potnuc_j,lroots_j,
     >                    nroots_j,iroot_j,mxroot,
     >                    nrs1_j,nrs2_j,nrs3_j,
     >                    nhole1_j,nelec3_j,ipt2_j,weight_j)

      if(debug)then
        write(6,*)' Information read from jobiph :'
        write(6,*)' ------------------------------'
        write(6,*)' nactel :',nactel_j
        write(6,*)' ispin  :',ispin_j
        write(6,*)' nsym   :',nsym_j
        write(6,*)' lsym   :',lSym_j
        write(6,*)' nfro   :',nfro_j
        write(6,*)' nish   :',nish_j
        write(6,*)' nash   :',nash_j
        write(6,*)' ndel   :',ndel_j
        write(6,*)' nbas   :',nbas_j
        write(6,*)' name   :'
        do 100 ii=1,mxorb
        if(name_j(ii)(4:4).eq.' ')write(6,*)name_j(ii)
100     continue
        write(6,*)' nconf  :',nconf_j
        write(6,*)' header :',header_j
        write(6,*)' title  :'
        do 200 ii=1,mxtit
        if(len_trim_cvb(title_j(ii)).gt.0)write(6,*)title_j(ii)
200     continue
        write(6,*)' potnuc :',potnuc_j
        write(6,*)' lroots :',lroots_j
        write(6,*)' nroots :',nroots_j
        write(6,*)' iroot  :',iroot_j
        write(6,*)' nrs1   :',nrs1_j
        write(6,*)' nrs2   :',nrs2_j
        write(6,*)' nrs3   :',nrs3_j
        write(6,*)' nhole1 :',nhole1_j
        write(6,*)' nelec3 :',nelec3_j
        write(6,*)' ipt2   :',ipt2_j
        write(6,*)' weight :',weight_j
      endif

      call daclos_cvb(lujob)
      Return
      End
