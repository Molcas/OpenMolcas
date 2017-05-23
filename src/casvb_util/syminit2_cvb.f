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
      subroutine syminit2_cvb(symelm,iorbrel,north,corth,
     >  irels,relorb,io,iorder,iorbs,a,b,
     >  rr,ri,vr,vi,intger,
     >  ifxorb,ifxstr,idelstr,
     >  iorts,irots,izeta)
      implicit real*8 (a-h,o-z)
      logical found
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension symelm(norb,norb,nsyme),iorbrel(ndimrel)
      dimension irels(2,*),relorb(norb,norb,*)
      dimension north(norb),corth(norb,*)
      dimension io(4,norbrel),iorder(norb,norbrel),iorbs(norb)
      dimension a(norb,norb),b(norb,norb)
      dimension rr(norb),ri(norb),vr(norb,norb),vi(norb,norb),
     >  intger(norb)
      dimension ifxorb(norb),ifxstr(nfxvb),idelstr(nzrvb)
      dimension iorts(2,nort),irots(2,ndrot),izeta(nsyme)
      save thresh
      data thresh/1.d-8/

c  Restore arrays :
      call rdioff_cvb(9,recinp,ioffs)
      call rdis_cvb(iorbrel,ndimrel,recinp,ioffs)
      call rdis_cvb(ifxorb,norb,recinp,ioffs)
      call rdis_cvb(ifxstr,nfxvb,recinp,ioffs)
      call rdis_cvb(idelstr,nzrvb,recinp,ioffs)
      call rdis_cvb(iorts,2*nort,recinp,ioffs)
      call rdis_cvb(irots,2*ndrot,recinp,ioffs)
      call rdis_cvb(izeta,nsyme,recinp,ioffs)

c  First check that minimum number of orbital relations has been given
c  (no cycle should be complete) :
      ishift=0
      do 10 ijrel=1,nijrel
      iorb=abs(iorbrel(1+ishift))
      jorb=abs(iorbrel(2+ishift))
      if(iorb.eq.jorb)goto 10
      call izero(iorbs,norb)
      iorbs(iorb)=1
      iorbs(jorb)=1
      ishift2=0
      do 30 ijrel2=1,nijrel
      iorb2=abs(iorbrel(1+ishift2))
      jorb2=abs(iorbrel(2+ishift2))
      if(ijrel2.eq.ijrel.or.iorb2.eq.jorb2)goto 30
      if(iorbs(iorb2).eq.1)then
        if(iorbs(jorb2).eq.1)then
          ncount=0
          do 40 i=1,norb
          if(iorbs(i).eq.1)then
            ncount=ncount+1
            intger(ncount)=i
          endif
40        continue
          write(6,'(2a,/,20i4)')' Too many orbital relations ',
     >      'involving orbitals :',(intger(ii),ii=1,ncount)
          write(6,'(a)')' Please reduce number of ORBREL cards.'
          call abend_cvb()
        else
          iorbs(jorb2)=1
        endif
      elseif(iorbs(jorb2).eq.1)then
        iorbs(iorb2)=1
      endif
30    ishift2=ishift2+3+iorbrel(3+ishift2)
10    ishift=ishift+3+iorbrel(3+ishift)

c  Diagonal orbital relations :
      nciorth=0
      call izero(north,norb)
      do 100 iorb=1,norb
c  Orbital conditions on IORB :
      call span0_cvb(norb,norb)
      ishift=0
      do 150 i=1,norbrel
      ior=iorbrel(1+ishift)
      jor=iorbrel(2+ishift)
      nrel=iorbrel(3+ishift)
      if(iorb.eq.ior.and.iorb.eq.jor)then
        call mxunit_cvb(b,norb)
        do 200 ir=nrel,1,-1
        irel=iorbrel(ir+3+ishift)
        call mxatb_cvb(symelm(1,1,irel),b,norb,norb,norb,a)
200     call fmove(a,b,norb*norb)
c  Everything that hasn't got eigenvalue +1 will be orthogonalised away
c  Unsymmetric diagonalisation :
        ifail=0
        call f02agf(b,norb,norb,rr,ri,vr,norb,vi,norb,intger,ifail)
        if(ifail.ne.0)then
          write(6,*)' Error in diagonalisation, IFAIL :',ifail
          call abend_cvb()
        endif
        do 250 ieig=1,norb
        if(abs(rr(ieig)-one).gt.thresh.or.abs(ri(ieig)).gt.thresh)then
          call addvec(vr(1,ieig),vr(1,ieig),vi(1,ieig),norb)
          call span1_cvb(vr(1,ieig),1,dum,norb,0)
        endif
250     continue
      endif
      ishift=ishift+3+nrel
150   continue
      if(plc_const)call rconstr_plc(iorb)
      call span2_cvb(corth(1,1+nciorth),north(iorb),dum,norb,0)
      nciorth=nciorth+north(iorb)
100   continue
      niorth=nciorth

c  Off-diagonal relations :
      call izero(io,2*norbrel)
      call izero(iorder,norb*norbrel)
      ijrel=0
      ishift=0
      do 300 i=1,norbrel
      iorb=iorbrel(1+ishift)
      jorb=iorbrel(2+ishift)
      nrel=iorbrel(3+ishift)
      if(iorb.ne.jorb)then
        ijrel=ijrel+1
        io(1,ijrel)=iorb
        io(2,ijrel)=jorb
        io(3,ijrel)=nrel
        io(4,ijrel)=3+ishift
        iorder(1,ijrel)=iorb
        iorder(2,ijrel)=jorb
        irel=iorbrel(4+ishift)
      endif
300   ishift=ishift+3+nrel
      nijrel=ijrel

c  Orbitals with constraints should be generating orbitals if possible:
      icnt=0
      do 400 iorb=1,norb
      if(north(iorb).gt.0)then
        icnt=icnt+1
        iorbs(icnt)=iorb
      endif
400   continue
      do 500 iorb=1,norb
      if(north(iorb).eq.0)then
        icnt=icnt+1
        iorbs(icnt)=iorb
      endif
500   continue
600   continue
c  Sort relations and define generating orbitals:
      do 700 i=1,norb
      iorb=iorbs(i)
      do 700 ii=1,nijrel
      if(iorder(1,ii).eq.iorb.or.iorder(2,ii).eq.iorb)then
c  Has IORB already been generated from KORB ? :
        do 800 j=1,i-1
        korb=iorbs(j)
800     if(iorder(1,ii).eq.korb.or.iorder(2,ii).eq.korb)goto 700
        if(iorder(1,ii).eq.iorb)then
          iorder(1,ii)=iorder(2,ii)
          iorder(2,ii)=iorb
        endif
        jorb=iorder(1,ii)
c  JORB will be generated from IORB
        if(north(jorb).ne.0)then
          write(6,'(2(a,i4),a)')' Attempting to generate orbital',jorb,
     >      ' from orbital',iorb,'  ---'
          write(6,'(2a,i4,a)')' the orbital conditions',
     >      ' for orbital',jorb,' cannot be enforced.'
          write(6,'(a)')' Please reduce number of ORBREL cards.'
          call abend_cvb()
        endif
        found=.false.
        do 900 jj=1,nijrel
        if(jj.ne.ii.and.(iorder(1,jj).eq.jorb.or.
     >    iorder(2,jj).eq.jorb))then
c  Is JORB involved in any other orbital relations ? :
          do 1000 j=1,i-1
          korb=iorbs(j)
1000      if(iorder(1,jj).eq.korb.or.iorder(2,jj).eq.korb)goto 900
          found=.true.
          if(iorder(1,jj).eq.jorb)then
            iorder(1,jj)=iorder(2,jj)
            iorder(2,jj)=jorb
          endif
          iorder(3,jj)=iorder(2,jj)
          call imove_cvb(iorder(3,ii),iorder(4,jj),norb-3)
c  KORB will be generated from IORB (via JORB) :
          iorder(2,jj)=iorder(2,ii)
        endif
900     continue
        if(found)goto 600
      endif
700   continue
c  Generate transformation matrix for each relation :
      do i=1,nijrel
      irels(1,i)=iorder(1,i)
      irels(2,i)=iorder(2,i)
      enddo
      do 1100 ijrel=1,nijrel
      call mxunit_cvb(relorb(1,1,ijrel),norb)
      do 1200 i=1,norb-1
      il=norb+2-i
      if(i.eq.1)il=2
      if(i.eq.norb)il=1
      if(iorder(il,ijrel).eq.0)goto 1200
      iorb=iorder(il,ijrel)
      j=i
1300  j=j+1
      jl=norb+2-j
      if(j.eq.1)jl=2
      if(j.eq.norb)jl=1
      if(iorder(jl,ijrel).eq.0)goto 1300
      jorb=iorder(jl,ijrel)
      do 1400 ii=1,nijrel
      if((iorb.eq.io(1,ii).or.iorb.eq.io(2,ii)).and.
     >   (jorb.eq.io(1,ii).or.jorb.eq.io(2,ii)))then
        nrel=io(3,ii)
        iaddr=io(4,ii)
c  Operate right-to-left :
        do 1500 ir=nrel,1,-1
        irel=iorbrel(ir+iaddr)
        if(jorb.eq.io(1,ii))then
          call mxatb_cvb(symelm(1,1,irel),relorb(1,1,ijrel),
     >      norb,norb,norb,a)
        else
          call mxattb_cvb(symelm(1,1,irel),relorb(1,1,ijrel),
     >      norb,norb,norb,a)
        endif
1500    call fmove(a,relorb(1,1,ijrel),norb*norb)
      endif
1400  continue
1200  continue
1100  continue
      return
      end
c  ********************
c  ** Symmetrization **
c  ********************
