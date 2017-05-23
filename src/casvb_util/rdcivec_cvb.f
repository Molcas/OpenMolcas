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
      Subroutine rdcivec_cvb(detvec,fn,reord)
************************************************************************
*                                                                      *
*     Read the contents of the JOBIPH file.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "rasdim.fh"
#include "jobiph_j.fh"
      character*(*) fn
      logical debug
      data debug/.false./
      dimension detvec(*)
      dimension ncix(8)
      logical reord

      iwr=0

      call getnci_cvb(ncix,nactel_j,ispin_j-1,lsym_j)
      ndet_j=ncix(1)

      lujob=15
      call daname_cvb(lujob,fn)
c Allocate at least NDET words for each vector, since this is
c required by csdtvc :
c      Call GetMem('OCIvec','Allo','Real',ipCI,nConf_j*nroots_j)
      Call GetMem('OCIvec','Allo','Real',ipCI,
     >  nConf_j*nroots_j+ndet_j-nconf_j)
      if(iwr.eq.0)then
        Do 200 i=1,nroots_j
        j=iroot_j(i)
        iDisk=iadr15_j(4)
        Do 300 k=1,j-1
300     Call dDaFile(LuJob,0,rdum,nConf_j,iDisk)
200     Call dDaFile(LuJob,2,Work(ipCI+(i-1)*nconf_j),
     >    nConf_j,iDisk)

        if(reord)then
          Call GetMem('ipci2','Allo','Real',ipCI2,nConf_j)
          call reord2_cvb(work(ipci),work(ipci2),1)
          call fmove(work(ipci2),work(ipci),nconf_j)
          Call GetMem('ipci2','Free','Real',ipCI2,idum)
        endif

        call csf2det_cvb(work(ipci),detvec,lsym_j,1)
      elseif(iwr.eq.1)then
        call csf2det_cvb(work(ipci),detvec,lsym_j,2)

        if(reord)then
          Call GetMem('ipci2','Allo','Real',ipCI2,nConf_j)
          call reord2_cvb(work(ipci),work(ipci2),0)
          call fmove(work(ipci2),work(ipci),nconf_j)
          Call GetMem('ipci2','Free','Real',ipCI2,idum)
        endif

        Do 400 i=1,nroots_j
        j=iroot_j(i)
        iDisk=iadr15_j(4)
        Do 500 k=1,j-1
500     Call dDaFile(LuJob,0,rdum,nConf_j,iDisk)
400     Call dDaFile(LuJob,1,Work(ipCI+(i-1)*nconf_j),
     >    nConf_j,iDisk)
      endif
      if (debug) then
        do 600 i=0,nroots_j-1
        write(6,'(a,i3,a)')' (CSF) CI vector ',i+1,' :'
        write(6,'(a)')' ---------------------'
600     call mxprint_cvb(work(ipci+nconf_j*i),1,nconf_j,0)
      endif
      Call GetMem('OCIvec','Free','Real',ipCI,idum)
      call daclos_cvb(lujob)
      Return
      end
      subroutine wrcivec_cvb(detvec,fn,reord)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "rasdim.fh"
#include "jobiph_j.fh"
      character*(*) fn
      logical debug
      data debug/.false./
      dimension detvec(*)
      dimension ncix(8)
      logical reord
      iwr=1

      call getnci_cvb(ncix,nactel_j,ispin_j-1,lsym_j)
      ndet_j=ncix(1)

      lujob=15
      call daname_cvb(lujob,fn)
c Allocate at least NDET words for each vector, since this is
c required by csdtvc :
c      Call GetMem('OCIvec','Allo','Real',ipCI,nConf_j*nroots_j)
      Call GetMem('OCIvec','Allo','Real',ipCI,
     >  nConf_j*nroots_j+ndet_j-nconf_j)
      if(iwr.eq.0)then
        Do 200 i=1,nroots_j
        j=iroot_j(i)
        iDisk=iadr15_j(4)
        Do 300 k=1,j-1
300     Call dDaFile(LuJob,0,rdum,nConf_j,iDisk)
200     Call dDaFile(LuJob,2,Work(ipCI+(i-1)*nconf_j),
     >    nConf_j,iDisk)

        if(reord)then
          Call GetMem('ipci2','Allo','Real',ipCI2,nConf_j)
          call reord2_cvb(work(ipci),work(ipci2),1)
          call fmove(work(ipci2),work(ipci),nconf_j)
          Call GetMem('ipci2','Free','Real',ipCI2,idum)
        endif

        call csf2det_cvb(work(ipci),detvec,lsym_j,1)
      elseif(iwr.eq.1)then
        call csf2det_cvb(work(ipci),detvec,lsym_j,2)

        if(reord)then
          Call GetMem('ipci2','Allo','Real',ipCI2,nConf_j)
          call reord2_cvb(work(ipci),work(ipci2),0)
          call fmove(work(ipci2),work(ipci),nconf_j)
          Call GetMem('ipci2','Free','Real',ipCI2,idum)
        endif

        Do 400 i=1,nroots_j
        j=iroot_j(i)
        iDisk=iadr15_j(4)
        Do 500 k=1,j-1
500     Call dDaFile(LuJob,0,rdum,nConf_j,iDisk)
400     Call dDaFile(LuJob,1,Work(ipCI+(i-1)*nconf_j),
     >    nConf_j,iDisk)
      endif
      if (debug) then
        do 600 i=0,nroots_j-1
        write(6,'(a,i3,a)')' (CSF) CI vector ',i+1,' :'
        write(6,'(a)')' ---------------------'
600     call mxprint_cvb(work(ipci+nconf_j*i),1,nconf_j,0)
      endif
      Call GetMem('OCIvec','Free','Real',ipCI,idum)
      call daclos_cvb(lujob)
      Return
      end
