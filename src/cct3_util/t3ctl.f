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
c     this file contains following routines:
c     nit3 (program)
c     getint
c        NB getint has been renamed to cct3_getint
c        to solve the conflict with lucia's getint
c     GetIntPoss
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine CCT3(ireturn)
c
c     this program calculate noniterative T3 contributions
c     to CCSD energy
c
#include "t31.fh"
#include "t32.fh"
c
c     work file declaration
       integer wrksize
#include "WrkSpc.fh"
c
c     variables
c
       real*8 eccsd,eaaa(1),eaab(1),eabb(1),ebbb(1),ec
       integer i,j,k,i3,j3,k3,symi,symj,symk,symi3,symj3,symk3
       integer symijk,symij
       integer nsg,keyyes,rc1,jup,ilow,posst,ssl1,ssl2,ssm1,ssm2
       integer symjstart,symjstop
       integer jstart,jstop,istart,istop
cpar
       integer id,counter
c
c       integer nhelp
c
#include "para_info.fh"
c
cpar
cstare Call SetTim
cstare call MPI_COMM_RANK(MPI_COMM_WORLD,myRank,rc)
cstare call MPI_COMM_SIZE(MPI_COMM_WORLD,nProcs,rc)

cI    **************  start   section **************
c
cI.1  read input data form INPDAT (reorg) and input file
      fullprint=0
      If (iPrintLevel(-1).LE.0) fullprint=-1
      If (fullprint.GE.0) then
         Write(6,*)
         Write(6,*)' **********************************'
         Write(6,*)'  Triples Contribution Calculation '
         Write(6,*)' **********************************'
         Write(6,*)
      EndIf
      call t3reainput
c
cI.2  write head to output file
      If (fullprint.GE.0) call t3wrhead
c
cI.3.1 calc. work space requirements to fix and help mediates
      call t3initfiles (wrksize)
      If (fullprint.GE.0)
     & write(6,*) ' Work space requirements :',wrksize
c
cI.3.2 allocate work space
c
      Call GetMem('CCT3','Max','Real',maxspace,maxspace)
      maxspace=maxspace-4
      if (maxspace.lt.wrksize) then
         write(6,*) ' Allocation of work space failed',maxspace
         write(6,*) ' Increase the value of the variable MOLCAS_MEM'
         Call Abend()
      end if
      Call GetMem('CCT3','Allo','Real',iOff,wrksize)
c
cI.3.3 set wrk = 0
      call cct3_mv0zero (wrksize,wrksize,Work(iOff))
      If (fullprint.GE.0) write(6,*) ' Allocation of work space : Done'
c
cI.4  read static integrals from INTSTA (reorg) file
      call t3reaintsta (Work(iOff),wrksize)
c
cI.5  divide fok to faa,fai,fii and dp
      call cct3_divfok (Work(iOff),wrksize,
     & mapdn,mapin,possn0,mapdp,mapip,possp0,
     & mapdfk1,mapifk1,possfk10,mapdfk2,mapifk2,possfk20,
     & mapdfk3,mapifk3,possfk30,mapdfk4,mapifk4,possfk40,
     & mapdfk5,mapifk5,possfk50,mapdfk6,mapifk6,possfk60,
     & mapddp1,mapidp1,possdp10,mapddp2,mapidp2,possdp20,rc1)
c
cI.6  read CCSD energy and amplitudes
      call t3reaccsd (Work(iOff),wrksize,
     & eccsd)
c
cI.7  get address vector T3IntPoss
c       they are located at the beggining of the t3nam file
      call GetIntPoss
c
cI.*  set energies=0
      eaaa=0.0d0
      eaab=0.0d0
      eabb=0.0d0
      ebbb=0.0d0
c
cI.par initialize paralell counter
      counter=0
c
c
c     ***** work section *****
c
c
      if (ijsegkey.eq.0) then
cNo Segmentation
         symimin=1
         symimax=nsym
      end if
c
cSegmented
cnoseg do 1400 symi=1,nsym
      do 1400 symi=symimin,symimax
       if (noa(symi).eq.0) then
       goto 1400
       end if
c
c      def symjstart,symjstop
c
       if (ijsegkey.eq.0) then
cNo Segmentation
          symjstart=1
          symjstop=symi
       else
c
cSegmented
       if (symimin.eq.symjmax) then
c      case symimin=symi=symimax
         symjstart=symjmin
         symjstop=symjmax
       else
c      case symimin<symjmax
c
         if(symi.eq.symimin) then
           symjstart=symjmin
           symjstop=symi
         else if (symi.eq.symimax) then
           symjstart=1
           symjstop=symjmax
         else
cc       sub case symimin<symi<symimax
           symjstart=1
           symjstop=symi
         end if
c
       end if
       end if
c
c
cnoseg do 1300 symj=1,symi
       do 1300 symj=symjstart,symjstop
       if ((symi.eq.symj).and.(noa(symi).le.1)) then
       goto 1300
       else if (noa(symj).eq.0) then
       goto 1300
       end if
c
c*    define sym(ij)
       symij=mmul(symi,symj)
c
c*    write symmetry starus of I and J
       if (fullprint.gt.0) then
       write(6,*) ' SYMI',symi,'SYMJ',symj
       end if
c
c
       if (symi.eq.symj) then
       ilow=2
       else
       ilow=1
       end if
c
c
c*    loop over inex I
c
c       define istart,istop
c
       if (ijsegkey.eq.0) then
cNo Segmentation
          istart=ilow
          istop=noa(symi)
       else
c
cSegmented
        if ((symi.eq.symimin).and.(symj.eq.symjmin)) then
          istart=imin
        else
          istart=ilow
        end if
c
        if ((symi.eq.symimax).and.(symj.eq.symjmax)) then
          istop=imax
        else
          istop=noa(symi)
        end if
c
       end if
c
cnoseg    do 1200 i=ilow,noa(symi)
          do 1200 i=istart,istop
c
c*    get integrals <ab|ic> for given i into R1(a,bc)
       call cct3_getint (Work(iOff),wrksize,
     & i,symi,possr10,mapdr1,mapir1,rc1)
c
c*    def upper limit for index j
       if (symi.eq.symj) then
       jup=i-1
       else
       jup=noa(symj)
       end if
c
c
c*    loop over index J
c
c*      def jstart,jstop
c
       if (ijsegkey.eq.0) then
cNo Segmentation
          jstart=1
          jstop=jup
       else
c
cSegmented
         if ((symimin.eq.symimax).and.(symjmin.eq.symjmax)) then
c        we are in the only symmetry block taken
c        into the consideration
           if (imin.eq.imax) then
              jstart=jmin
              jstop=jmax
           else if (i.eq.imin) then
              jstart=jmin
              jstop=jup
           else if (i.eq.imax) then
              jstart=1
              jstop=jmax
           else
              jstart=1
              jstop=jup
           end if
         else if ((symi.eq.symimin).and.(symj.eq.symjmin)) then
c        we are in initial symmetry block
           if (i.eq.imin) then
             jstart=jmin
             jstop=jup
           else
             jstart=1
             jstop=jup
           end if
         else if ((symi.eq.symimax).and.(symj.eq.symjmax)) then
c        we are in terminal symmetry block
           if (i.eq.imax) then
             jstart=1
             jstop=jmax
           else
             jstart=1
             jstop=jup
           end if
         else
c        we are in intermediate symmetry block
           jstart=1
           jstop=jup
         end if
c
       end if
c
c
cnoseg do 1200 j=1,jup
       do 1201 j=jstart,jstop
c
c*    get integrals <ab|jc> for given j into R2(a,bc)
       call cct3_getint (Work(iOff),wrksize,
     & j,symj,possr20,mapdr2,mapir2,rc1)
c
c
       do 1100 symk=1,nsym
       if (noa(symk).eq.0) then
       goto 1100
       end if
c
       if (fullprint.gt.1) then
       write(6,*) ' SYMI',symi,'SYMJ',symj,'SYMK',symk
       end if
c
c*    define sym(ijk)
       symijk=mmul(symij,symk)
c
c
c*    loop over index K
c
       do 1000 k=1,noa(symk)
       if (fullprint.ge.2) then
       write (6,999) i,j,k
999    format (' I,J,K ',3(i3,1x))
       end if
c
cpar    update paralell counter, choose proper id for this 'portion'
c       and skip if this portion is not for myRank
        counter=counter+1
        id=mod(counter,nProcs)
        if (myRank.ne.id) goto 1000
c
c*    get integrals <ab|kc> for given k into R3(a,bc)
       call cct3_getint (Work(iOff),wrksize,
     & k,symk,possr30,mapdr3,mapir3,rc1)
c
c
c1    ***** aaa spin combination *****
c
c
c1*   def keyyes
       if (symj.gt.symk) then
       keyyes=1
       else if (symj.eq.symk) then
       if (j.gt.k) then
       keyyes=1
       else
       keyyes=0
       end if
       else
       keyyes=0
       end if
c
       if (keyyes.eq.1) then
c
c1.*  define maps of W(abc)
       call cct3_t3grc0(3,5,3,3,3,0,symijk,possw0,posst,mapdw,mapiw)
c
c1.*  vanish W(abc)
       call stz (Work(iOff),wrksize,
     & mapdw)
c
c
c1.1  permutations (ijk) P(a,bc) (general sign+)
       nsg=1
c
c
c1.1.1V graph
c
c1.1.1*            def L1(bc,d) <- R3(b,cd) for given k
       call defv (Work(iOff),wrksize,
     & 1,possl10,mapdl1,mapil1,ssl1,
     & mapdr3,mapir3,symk,rc1)
c
c1.1.1*            ext M1(da) <- T2aaaa(da,ij) for given i,j
       call ext(Work(iOff),wrksize,
     & 4,7,i,j,0,symi,symj,0,mapdt21,mapit21,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c1.1.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,
     & rc1)
       ssm2=ssm1
c
c1.1.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c1.1.1*            pack W(abc) <-  P(a,bc) [L2(bc,a)] (minus due to the
c     usung Tikda instead of T2ikad)
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
c
c1.1.2O graph
c
c1.1.2*            ext L1(bc,l) <- T2aaaa(bc,kl) for given k
       call ext(Work(iOff),wrksize,
     & 4,3,k,0,0,symk,0,0,mapdt21,mapit21,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c1.1.2*     ext M1(l,a) <- W11(l,a,ij)=<la||ij>aaaa for given ij
       call ext(Work(iOff),wrksize,
     & 4,7,i,j,0,symi,symj,0,mapdw11,mapiw11,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c1.1.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c1.1.2*            pack W(abc) <- P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c
c1.2  permutations (ikj) P(a,bc)
       nsg=-1
c
c
c1.2.1V graph
c
c1.2.1*            def L1(bc,d) <- R2(b,cd) for given j
       call defv (Work(iOff),wrksize,
     & 1,possl10,mapdl1,mapil1,ssl1,
     & mapdr2,mapir2,symj,rc1)
c
c1.2.1*            ext M1(da) <- T2aaaa(da,ik) for given i,k
       call ext(Work(iOff),wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdt21,mapit21,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c1.2.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,
     & rc1)
       ssm2=ssm1
c
c1.2.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c1.2.1*            pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to
c     usung Tikda instead of T2ikad)
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
c
c1.2.2O graph
c
c1.2.2*            ext L1(bc,l) <- T2aaaa(bc,jl) for given j
       call ext(Work(iOff),wrksize,
     & 4,3,j,0,0,symj,0,0,mapdt21,mapit21,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c1.2.2*     ext M1(l,a) <- W11(l,a,ik)=<la||ik>aaaa for given ik
       call ext(Work(iOff),wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdw11,mapiw11,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c1.2.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c1.2.2*            pack W(abc) <- P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c
c1.3  permutations (jki) P(a,bc)
       nsg=1
c
c
c1.3.1V graph
c
c1.3.1*            def L1(bc,d) <- R1(b,cd) for given i
       call defv (Work(iOff),wrksize,
     & 1,possl10,mapdl1,mapil1,ssl1,
     & mapdr1,mapir1,symi,rc1)
c
c1.3.1*            ext M1(da) <- T2aaaa(da,jk) for given j,k
       call ext(Work(iOff),wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdt21,mapit21,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c1.3.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,
     & rc1)
       ssm2=ssm1
c
c1.3.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c1.3.1*            pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to
c     usung Tjkda instead of T2jkad)
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
c
c1.3.2O graph
c
c1.3.2*            ext L1(bc,l) <- T2aaaa(bc,kl) for given i
       call ext(Work(iOff),wrksize,
     & 4,3,i,0,0,symi,0,0,mapdt21,mapit21,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c1.3.2*     ext M1(l,a) <- W11(l,a,jk)=<la||jk>aaaa for given jk
       call ext(Work(iOff),wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdw11,mapiw11,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c1.3.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c1.3.2*            pack W(abc) <-  P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c
c1.4  add singles
c
c1.4.0mov V <- -W
       call cct3_t3grc0(3,5,3,3,3,0,symijk,possv0,posst,mapdv,mapiv)
       call minusa (Work(iOff),wrksize,
     & mapdw,-1.0d0)
       call setb (Work(iOff),wrksize,
     & mapdw,mapdv,1.0d0)
c
       if (typt3.gt.1) then
c1.4.1add part W2 . T1
       call t3sgl(Work(iOff),wrksize,
     & mapdv,symijk,mapdt11,mapit11,mapdt12,mapit12,
     & mapdw21,mapiw21,mapdw22,mapiw22,
     & 1,i,j,k,symi,symj,symk,rc1,
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
c
       if (typt3.eq.3) then
c1.4.2add part T2 . U
       call t3sgl(Work(iOff),wrksize,
     & mapdv,symijk,mapdfk3,mapifk3,mapdfk4,mapifk4,
     & mapdt21,mapit21,mapdt22,mapit22,
     & 1,i,j,k,symi,symj,symk,rc1,
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
c
c
c1.5  divide by denominators and calc energy contribution
c
c
c1.5.1divide by den.
       call t3div(Work(iOff),wrksize,
     & mapdw,mapdv,symijk,mapddp1,mapidp1,mapddp2,
     & mapidp2,1,i,j,k,symi,symj,symk,ec,rc1)
c
c1.5.2add energy contribution
       eaaa=eaaa+ec
c
       end if
c
c
c2    ***** aab spin combination *****
c
c
c2*   def keyyes
       if (k.le.nob(symk)) then
       keyyes=1
       else
       keyyes=0
       end if
c
       if (keyyes.eq.1) then
c
c2.*  define maps of W(abc)
       call cct3_t3grc0(3,1,3,3,4,0,symijk,possw0,posst,mapdw,mapiw)
c
c2.*  vanish W(abc)
       call stz (Work(iOff),wrksize,
     & mapdw)
c
c
c2.1  permutations (ijk) P(a,b) (c)
       nsg=1
c
c2.1.1V graph
c
c2.1.1*            def L1(b,c,d)aba <- R3(b,cd) for given k
       call defv (Work(iOff),wrksize,
     & 4,possl10,mapdl1,mapil1,ssl1,
     & mapdr3,mapir3,symk,rc1)
c
c2.1.1*            ext M1(da) <- T2aaaa(da,ij) for given i,j
       call ext(Work(iOff),wrksize,
     & 4,7,i,j,0,symi,symj,0,mapdt21,mapit21,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.1.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,
     & rc1)
       ssm2=ssm1
c
c2.1.1*            mult L2(b,c,a) <- L1(b,c,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.1.1*            pack W(ab,c) <- - P(a,b) [L2(b,c,a)] (minus is due to
c     usung Tijda instead of T2ijad)
       call t3addpck(Work(iOff),wrksize,
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
c
c2.1.2O graph
c
c2.1.2*            ext L1(b,c,l) <- T2abab(b,c,l,k) for given k
       call ext(Work(iOff),wrksize,
     & 4,4,k,0,0,symk,0,0,mapdt23,mapit23,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c2.1.2*     ext M1(l,a) <- W12(l,a,ij)=<la||ij>aaaa for given ij
       call ext(Work(iOff),wrksize,
     & 4,7,i,j,0,symi,symj,0,mapdw11,mapiw11,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.1.2*            mult L2(b,c,a) <- L1(b,c,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.1.2*            pack W(ab,c) <- - P(a,b) [L2(b,c,a)] (-, premutation in V)
       call t3addpck(Work(iOff),wrksize,
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
c
c
c2.2  permutations (ijk) (cab) do not nontribute
c
c
c2.3  permutations (ikj) P(a,b) (c)
       nsg=-1
c
c2.3.1V graph
c
c2.3.1*            def L1(b,c,d)abb <- R2(b,cd) for given j
       call defv (Work(iOff),wrksize,
     & 3,possl10,mapdl1,mapil1,ssl1,
     & mapdr2,mapir2,symj,rc1)
c
c2.3.1*            ext M1(a,d) <- T2abab(a,d,i,k) for given i,k
       call ext(Work(iOff),wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdt23,mapit23,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.3.1*            map M2(d,a) <- M1(a,d)
       call cct3_map (Work(iOff),wrksize,
     & 2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,
     & possm20,posst,rc1)
       ssm2=ssm1
c
c2.3.1*            mult L2(b,c,a) <- L1(b,c,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.3.1*            pack W(ab,c) <-  P(a,b) [L2(b,c,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c2.3.2O graph
c
c2.3.2*            ext L1(b,c,l) <- T2abab(b,c,j,l) for given j
       call ext(Work(iOff),wrksize,
     & 4,3,j,0,0,symj,0,0,mapdt23,mapit23,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c2.3.2*     ext M1(l,a) <- W14(l,a,ik)=<la||ik>baab for given ik
       call ext(Work(iOff),wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdw14,mapiw14,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.3.2*            mult L2(b,c,a) <- L1(b,c,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.3.2*            pack W(ab,c) <- P(a,b) [L2(b,c,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c
c2.4  permutations (ikj) (cab)
       nsg=-1
c
c
c2.4.1V graph
c
c2.4.1*            def L1(ab,d)aaa <- R2(a,bc) for given j
       call defv (Work(iOff),wrksize,
     & 1,possl10,mapdl1,mapil1,ssl1,
     & mapdr2,mapir2,symj,rc1)
c
c2.4.1*            ext M1(d,c) <- T2abab(d,c,i,k) for given i,k
       call ext(Work(iOff),wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdt23,mapit23,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.4.1*            mult L2(ab,c) <- L1(ab,d) . M1(d,c)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.4.1*            add (pack) W(ab,c) <-  [L2(b,c,a)] (- due to permuted T)
       call cct3_add (Work(iOff),wrksize,
     & 3,3,0,0,0,0,1,1,-1.0d0*nsg,mapdl2,symijk,
     & mapdw,mapiw,symijk,rc1)
c
c2.4.2O graph
c
c2.4.2*            ext L1(ab,l) <- T2aaaa(ab,jl) for given j
       call ext(Work(iOff),wrksize,
     & 4,3,j,0,0,symj,0,0,mapdt21,mapit21,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c2.4.2*     ext M1(l,c) <- W13(l,c,ik)=<lc||ik>abab for given ik
       call ext(Work(iOff),wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdw13,mapiw13,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.4.2*            mult L2(ab,c) <- L1(ab,l) . M1(l,c)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.4.2*            add (pack) W(ab,c) <- [L2(ab,c)]
       call cct3_add (Work(iOff),wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl2,symijk,
     & mapdw,mapiw,symijk,rc1)
c
c
c2.5  permutations (jki) P(a,b) (c)
       nsg=1
c
c
c2.5.1V graph
c
c2.5.1*            def L1(b,c,d)abb <- R1(b,cd) for given i
       call defv (Work(iOff),wrksize,
     & 3,possl10,mapdl1,mapil1,ssl1,
     & mapdr1,mapir1,symi,rc1)
c
c2.5.1*            ext M1(a,d) <- T2abab(a,d,j,k) for given j,k
       call ext(Work(iOff),wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdt23,mapit23,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.5.1*            map M2(d,a) <- M1(a,d)
       call cct3_map (Work(iOff),wrksize,
     & 2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,
     & possm20,posst,rc1)
       ssm2=ssm1
c
c2.5.1*            mult L2(b,c,a) <- L1(b,c,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.5.1*            pack W(ab,c) <-  P(a,b) [L2(b,c,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c2.5.2O graph
c
c2.5.2*            ext L1(b,c,l) <- T2abab(b,c,i,l) for given i
       call ext(Work(iOff),wrksize,
     & 4,3,i,0,0,symi,0,0,mapdt23,mapit23,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c2.5.2*     ext M1(l,a) <- W14(l,a,jk)=<la||jk>baab for given jk
       call ext(Work(iOff),wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdw14,mapiw14,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.5.2*            mult L2(b,c,a) <- L1(b,c,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.5.2*            pack W(ab,c) <- P(a,b) [L2(b,c,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c
c2.6  permutations (ikj) (cab)
       nsg=1
c
c
c2.6.1V graph
c
c2.6.1*            def L1(ab,d)aaa <- R1(a,bc) for given i
       call defv (Work(iOff),wrksize,
     & 1,possl10,mapdl1,mapil1,ssl1,
     & mapdr1,mapir1,symi,rc1)
c
c2.6.1*            ext M1(d,c) <- T2abab(d,c,j,k) for given j,k
       call ext(Work(iOff),wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdt23,mapit23,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.6.1*            mult L2(ab,c) <- L1(ab,d) . M1(d,c)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.6.1*            add (pack) W(ab,c) <-  [L2(b,c,a)] (- due to permuted T)
       call cct3_add (Work(iOff),wrksize,
     & 3,3,0,0,0,0,1,1,-1.0d0*nsg,mapdl2,symijk,
     & mapdw,mapiw,symijk,rc1)
c
c2.6.2O graph
c
c2.6.2*            ext L1(ab,l) <- T2aaaa(ab,il) for given i
       call ext(Work(iOff),wrksize,
     & 4,3,i,0,0,symi,0,0,mapdt21,mapit21,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c2.6.2*     ext M1(l,c) <- W13(l,c,jk)=<lc||jk>abab for given jk
       call ext(Work(iOff),wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdw13,mapiw13,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c2.6.2*            mult L2(ab,c) <- L1(ab,l) . M1(l,c)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c2.6.2*            add (pack) W(ab,c) <- [L2(ab,c)]
       call cct3_add (Work(iOff),wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl2,symijk,
     & mapdw,mapiw,symijk,rc1)
c
c
c2.7  add singles
c
c
c2.7.0mov V <- W
       call cct3_t3grc0(3,1,3,3,4,0,symijk,possv0,posst,mapdv,mapiv)
       call minusa (Work(iOff),wrksize,
     & mapdw,-1.0d0)
       call setb (Work(iOff),wrksize,
     & mapdw,mapdv,1.0d0)
c
       if (typt3.gt.1) then
c2.7.1add part W2 . T1
       call t3sgl(Work(iOff),wrksize,
     & mapdv,symijk,mapdt11,mapit11,mapdt12,mapit12,
     & mapdw21,mapiw21,mapdw23,mapiw23,
     & 2,i,j,k,symi,symj,symk,rc1,
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
c
       if (typt3.eq.3) then
c2.7.2add part T2 . U
       call t3sgl(Work(iOff),wrksize,
     & mapdv,symijk,mapdfk3,mapifk3,mapdfk4,mapifk4,
     & mapdt21,mapit21,mapdt23,mapit23,
     & 2,i,j,k,symi,symj,symk,rc1,
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
c
c
c2.8  divide by denominators and calc energy contribution
c
c
c2.8.1divide by den.
       call t3div(Work(iOff),wrksize,
     & mapdw,mapdv,symijk,mapddp1,mapidp1,mapddp2,
     & mapidp2,2,i,j,k,symi,symj,symk,ec,rc1)
c
c2.8.2add energy contribution
       eaab=eaab+ec
c
       end if
c
c
c3    ***** abb spin combination *****
c
c
c     Note:
c     in spin combination 3-abb indexes are changed as follows
c
       i3=k
       symi3=symk
       j3=i
       symj3=symi
       k3=j
       symk3=symj
c
c     therefore, also R files are mixed:
c     R13 = R3, R23=R1 and R33=R2, symijk3=symijk
c
c3.*  def keyyes
       if (symj3.gt.symk3) then
       if ((j3.le.nob(symj3)).and.(k3.le.nob(symk3))) then
       keyyes=1
       else
       keyyes=0
       end if
       else if (symj3.eq.symk3) then
       if ((j3.gt.k3).and.(j3.le.nob(symj3))) then
       keyyes=1
       else
       keyyes=0
       end if
       else
       keyyes=0
       end if
c
       if (keyyes.eq.1) then
c
c3.*  define maps of W(a,bc)
       call cct3_t3grc0(3,2,3,4,4,0,symijk,possw0,posst,mapdw,mapiw)
c
c3.*  vanish W(abc)
       call stz (Work(iOff),wrksize,
     & mapdw)
c
c
c3.1  permutation (ijk) (abc)
       nsg=1
c
c
c3.1.1V graph
c
c3.1.1*            def L1(bc,d) <- R33(b,cd) for given k3
       call defv (Work(iOff),wrksize,
     & 2,possl10,mapdl1,mapil1,ssl1,
     & mapdr2,mapir2,symk3,rc1)
c
c3.1.1*            ext M1(a,d) <- T2abab(a,d,i3,j3) for given i3,j3
       call ext(Work(iOff),wrksize,
     & 4,7,i3,j3,0,symi3,symj3,0,mapdt23,mapit23,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.1.1*            map M2(d,a) <- M1(a,d)
       call cct3_map (Work(iOff),wrksize,
     & 2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,
     & possm20,posst,rc1)
       ssm2=ssm1
c
c3.1.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.1.1map L1(a,bc) <- L2(bc,a)
       call cct3_map (Work(iOff),wrksize,
     & 3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,
     & possl10,posst,rc1)
       ssl1=ssl2
c3.1.1*            add (pack) W(a,bc) <- [L1(a,bc)]
       call cct3_add (Work(iOff),wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl1,symijk,
     & mapdw,mapiw,symijk,rc1)
c
c3.1.2O graph
c
c3.1.2*            ext L1(bc,l) <- T2bbbb(bc,k3l) for given k3
       call ext(Work(iOff),wrksize,
     & 4,3,k3,0,0,symk3,0,0,mapdt22,mapit22,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c3.1.2*     ext M1(l,a) <- W14(l,a,i3j3)=<la||i3j3>baab for given i3,j3
       call ext(Work(iOff),wrksize,
     & 4,7,i3,j3,0,symi3,symj3,0,mapdw14,mapiw14,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.1.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.1.2map L1(a,bc) <- L2(bc,a)
       call cct3_map (Work(iOff),wrksize,
     & 3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,
     & possl10,posst,rc1)
       ssl1=ssl2

c3.1.2*            add (pack) W(a,bc) <- [L1(a,bc)]
       call cct3_add (Work(iOff),wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl1,symijk,
     & mapdw,mapiw,symijk,rc1)
c
c
c3.2  permutations (ijk)(bac),(ijk)(cab)
       nsg=1
c
c
c3.2.1V graph
c
c3.2.1*            def L1(a,c,d)aba <- R33(a,cd) for given k3
       call defv (Work(iOff),wrksize,
     & 4,possl10,mapdl1,mapil1,ssl1,
     & mapdr2,mapir2,symk3,rc1)
c
c3.2.1*            ext M1(d,b) <- T2abab(d,b,i3,j3) for given i3,j3
       call ext(Work(iOff),wrksize,
     & 4,7,i3,j3,0,symi3,symj3,0,mapdt23,mapit23,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.2.1*            mult L2(ac,b) <- L1(a,c,d) . M1(d,b)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.2.1*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm T)
       call t3addpck (Work(iOff),wrksize,
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,
     & 0,rc1)
c
c3.2.2O graph
c
c3.2.2*            ext L1(a,c,l) <- T2abab(a,c,l,k3) for given k3
       call ext(Work(iOff),wrksize,
     & 4,4,k3,0,0,symk3,0,0,mapdt23,mapit23,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c3.2.2*     ext M1(l,b) <- W13(l,b,i3j3)=<la||i3j3>abab for given i3,j3
       call ext(Work(iOff),wrksize,
     & 4,7,i3,j3,0,symi3,symj3,0,mapdw13,mapiw13,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.2.2*            mult L2(ac,b) <- L1(a,c,l) . M1(l,b)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.2.2*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm V)
       call t3addpck (Work(iOff),wrksize,
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,
     & 0,rc1)
c
c
c3.3  permutation (ikj) (abc)
       nsg=-1
c
c
c3.3.1V graph
c
c3.3.1*            def L1(bc,d) <- R23(b,cd) for given j3
       call defv (Work(iOff),wrksize,
     & 2,possl10,mapdl1,mapil1,ssl1,
     & mapdr1,mapir1,symj3,rc1)
c
c3.3.1*            ext M1(a,d) <- T2abab(a,d,i3,k3) for given i3,k3
       call ext(Work(iOff),wrksize,
     & 4,7,i3,k3,0,symi3,symk3,0,mapdt23,mapit23,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.3.1*            map M2(d,a) <- M1(a,d)
       call cct3_map (Work(iOff),wrksize,
     & 2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,
     & possm20,posst,rc1)
       ssm2=ssm1
c
c3.3.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.3.1map L1(a,bc) <- L2(bc,a)
       call cct3_map (Work(iOff),wrksize,
     & 3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,
     & possl10,posst,rc1)
       ssl1=ssl2

c3.3.1*            add (pack) W(a,bc) <- [L1(a,bc)]
       call cct3_add (Work(iOff),wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl1,symijk,
     & mapdw,mapiw,symijk,rc1)
c
c3.3.2O graph
c
c3.3.2*            ext L1(bc,l) <- T2bbbb(bc,j3l) for given j3
       call ext(Work(iOff),wrksize,
     & 4,3,j3,0,0,symj3,0,0,mapdt22,mapit22,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c3.3.2*     ext M1(l,a) <- W14(l,a,i3k3)=<la||i3k3>baab for given i3,k3
       call ext(Work(iOff),wrksize,
     & 4,7,i3,k3,0,symi3,symk3,0,mapdw14,mapiw14,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.3.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.3.2map L1(a,bc) <- L2(bc,a)
       call cct3_map (Work(iOff),wrksize,
     & 3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,
     & possl10,posst,rc1)
       ssl1=ssl2

c3.3.2*            add (pack) W(a,bc) <- [L1(a,bc)]
       call cct3_add (Work(iOff),wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl1,symijk,
     & mapdw,mapiw,symijk,rc1)
c
c
c3.4  permutations (ikj)(bac),(ikj)(cab)
       nsg=-1
c
c
c3.4.1V graph
c
c3.4.1*            def L1(a,c,d)aba <- R23(a,cd) for given j3
       call defv (Work(iOff),wrksize,
     & 4,possl10,mapdl1,mapil1,ssl1,
     & mapdr1,mapir1,symj3,rc1)
c
c3.4.1*            ext M1(d,b) <- T2abab(d,b,i3,k3) for given i3,k3
       call ext(Work(iOff),wrksize,
     & 4,7,i3,k3,0,symi3,symk3,0,mapdt23,mapit23,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.4.1*            mult L2(ac,b) <- L1(a,c,d) . M1(d,b)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.4.1*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm T)
       call t3addpck (Work(iOff),wrksize,
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,
     & 0,rc1)
c
c3.4.2O graph
c
c3.4.2*            ext L1(a,c,l) <- T2abab(a,c,l,j3) for given j3
       call ext(Work(iOff),wrksize,
     & 4,4,j3,0,0,symj3,0,0,mapdt23,mapit23,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c3.4.2*     ext M1(l,b) <- W13(l,b,i3k3)=<la||i3k3>abab for given i3,k3
       call ext(Work(iOff),wrksize,
     & 4,7,i3,k3,0,symi3,symk3,0,mapdw13,mapiw13,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.4.2*            mult L2(ac,b) <- L1(a,c,l) . M1(l,b)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.4.2*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm V)
       call t3addpck (Work(iOff),wrksize,
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,
     & 0,rc1)
c
c
c3.5  permutations (jki)(abc) do not contribute
c
c
c3.6  permutations (jki)(bac),(jki)(cab)
       nsg=1
c
c
c3.6.1V graph
c
c3.6.1*            def L1(a,c,d)abb <- R13(a,cd) for given i3
       call defv (Work(iOff),wrksize,
     & 3,possl10,mapdl1,mapil1,ssl1,
     & mapdr3,mapir3,symi3,rc1)
c
c3.6.1*            ext M1(db) <- T2bbbb(db,j3k3) for given j3,k3
       call ext(Work(iOff),wrksize,
     & 4,7,j3,k3,0,symj3,symk3,0,mapdt22,mapit22,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.6.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,
     & rc1)
       ssm2=ssm1
c
c3.6.1*            mult L2(ac,b) <- L1(a,c,d) . M2(d,b)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.6.1*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm T)
       call t3addpck (Work(iOff),wrksize,
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,
     & 0,rc1)
c
c3.6.2O graph
c
c3.6.2*            ext L1(a,c,l) <- T2abab(a,c,i3,l) for given i3
       call ext(Work(iOff),wrksize,
     & 4,3,i3,0,0,symi3,0,0,mapdt23,mapit23,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c3.6.2*     ext M1(l,b) <- W12(l,b,j3k3)=<la||j3k3>bbbb for given j3,k3
       call ext(Work(iOff),wrksize,
     & 4,7,j3,k3,0,symj3,symk3,0,mapdw12,mapiw12,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c3.6.2*            mult L2(ac,b) <- L1(a,c,l) . M1(l,b)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c3.6.2*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)]
       call t3addpck (Work(iOff),wrksize,
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,
     & 0,rc1)
c
c
c3.7  add singles
c
c
c3.7.0mov V <- W
       call cct3_t3grc0(3,2,3,4,4,0,symijk,possv0,posst,mapdv,mapiv)
       call minusa (Work(iOff),wrksize,
     & mapdw,-1.0d0)
       call setb (Work(iOff),wrksize,
     & mapdw,mapdv,1.0d0)
c
       if (typt3.gt.1) then
c3.7.1add part W2 . T1
       call t3sgl(Work(iOff),wrksize,
     & mapdv,symijk,mapdt11,mapit11,mapdt12,mapit12,
     & mapdw23,mapiw23,mapdw22,mapiw22,
     & 3,i3,j3,k3,symi3,symj3,symk3,rc1,
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
c
       if (typt3.eq.3) then
c3.7.2add part T2 . U
       call t3sgl(Work(iOff),wrksize,
     & mapdv,symijk,mapdfk3,mapifk3,mapdfk4,mapifk4,
     & mapdt23,mapit23,mapdt22,mapit22,
     & 3,i3,j3,k3,symi3,symj3,symk3,rc1,
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
c
c
c3.8  divide by denominators and calc energy contribution
c
c
c3.8.1divide by den.
       call t3div(Work(iOff),wrksize,
     & mapdw,mapdv,symijk,mapddp1,mapidp1,mapddp2,
     & mapidp2,3,i3,j3,k3,symi3,symj3,symk3,ec,rc1)
c
c3.8.2add energy contribution
       eabb=eabb+ec
c
       end if
c
c
c4    ***** bbb spin combination *****
c
c
c4.*  def keyyes
       if (symj.gt.symk) then
       keyyes=1
       else if (symj.eq.symk) then
       if (j.gt.k) then
       keyyes=1
       else
       keyyes=0
       end if
       else
       keyyes=0
       end if
c
       if ((i.gt.nob(symi)).or.(j.gt.nob(symj))
     & .or.(k.gt.nob(symk))) then
       keyyes=0
       end if
c
       if (keyyes.eq.1) then
c
c4.*  define maps of W(abc)
       call cct3_t3grc0(3,5,4,4,4,0,symijk,possw0,posst,mapdw,mapiw)
c
c4.*  vanish W(abc)
       call stz (Work(iOff),wrksize,
     & mapdw)
c
c
c4.1  permutations (ijk) P(a,bc) (general sign+)
       nsg=1
c
c
c4.1.1V graph
c
c4.1.1*            def L1(bc,d) <- R3(b,cd) for given k
       call defv (Work(iOff),wrksize,
     & 2,possl10,mapdl1,mapil1,ssl1,
     & mapdr3,mapir3,symk,rc1)
c
c4.1.1*            ext M1(da) <- T2bbbb(da,ij) for given i,j
       call ext(Work(iOff),wrksize,
     & 4,7,i,j,0,symi,symj,0,mapdt22,mapit22,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c4.1.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,
     & rc1)
       ssm2=ssm1
c
c4.1.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c4.1.1*            pack W(abc) <-  P(a,bc) [L2(bc,a)] (minus due to the
c     usung Tikda instead of T2ikad)
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
c
c4.1.2O graph
c
c4.1.2*            ext L1(bc,l) <- T2bbbb(bc,kl) for given k
       call ext(Work(iOff),wrksize,
     & 4,3,k,0,0,symk,0,0,mapdt22,mapit22,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c4.1.2*     ext M1(l,a) <- W12(l,a,ij)=<la||ij>bbbb for given ij
       call ext(Work(iOff),wrksize,
     & 4,7,i,j,0,symi,symj,0,mapdw12,mapiw12,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c4.1.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c4.1.2*            pack W(abc) <- P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c
c4.2  permutations (ikj) P(a,bc)
       nsg=-1
c
c
c4.2.1V graph
c
c4.2.1*            def L1(bc,d) <- R2(b,cd) for given j
       call defv (Work(iOff),wrksize,
     & 2,possl10,mapdl1,mapil1,ssl1,
     & mapdr2,mapir2,symj,rc1)
c
c4.2.1*            ext M1(da) <- T2bbbb(da,ik) for given i,k
       call ext(Work(iOff),wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdt22,mapit22,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c4.2.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,
     & rc1)
       ssm2=ssm1
c
c4.2.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c4.2.1*            pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to
c     usung Tikda instead of T2ikad)
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
c
c4.2.2O graph
c
c4.2.2*            ext L1(bc,l) <- T2bbbb(bc,jl) for given j
       call ext(Work(iOff),wrksize,
     & 4,3,j,0,0,symj,0,0,mapdt22,mapit22,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c4.2.2*     ext M1(l,a) <- W12(l,a,ik)=<la||ik>bbbb for given ik
       call ext(Work(iOff),wrksize,
     & 4,7,i,k,0,symi,symk,0,mapdw12,mapiw12,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c4.2.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c4.2.2*            pack W(abc) <- P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c
c4.3  permutations (jki) P(a,bc)
       nsg=1
c
c
c4.3.1V graph
c
c4.3.1*            def L1(bc,d) <- R1(b,cd) for given i
       call defv (Work(iOff),wrksize,
     & 2,possl10,mapdl1,mapil1,ssl1,
     & mapdr1,mapir1,symi,rc1)
c
c4.3.1*            ext M1(da) <- T2bbbb(da,jk) for given j,k
       call ext(Work(iOff),wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdt22,mapit22,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c4.3.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,
     & rc1)
       ssm2=ssm1
c
c4.3.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c4.3.1*            pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to
c     usung Tjkda instead of T2jkad)
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
c
c4.3.2O graph
c
c4.3.2*            ext L1(bc,l) <- T2bbbb(bc,kl) for given i
       call ext(Work(iOff),wrksize,
     & 4,3,i,0,0,symi,0,0,mapdt22,mapit22,1,
     & possl10,mapdl1,mapil1,ssl1,rc1)
c
c4.3.2*     ext M1(l,a) <- W12(l,a,jk)=<la||jk>bbbb for given jk
       call ext(Work(iOff),wrksize,
     & 4,7,j,k,0,symj,symk,0,mapdw12,mapiw12,1,
     & possm10,mapdm1,mapim1,ssm1,rc1)
c
c4.3.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,
     & mapdl2,mapil2,ssl2,possl20,rc1)
c
c4.3.2*            pack W(abc) <-  P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
c
c
c4.4  add singles
c
c
c4.4.0mov V <- W
       call cct3_t3grc0(3,5,4,4,4,0,symijk,possv0,posst,mapdv,mapiv)
       call minusa (Work(iOff),wrksize,
     & mapdw,-1.0d0)
       call setb (Work(iOff),wrksize,
     & mapdw,mapdv,1.0d0)
c
       if (typt3.gt.1) then
c4.4.1add part W2 . T1
       call t3sgl(Work(iOff),wrksize,
     & mapdv,symijk,mapdt12,mapit12,mapdt11,mapit11,
     & mapdw22,mapiw22,mapdw21,mapiw21,
     & 1,i,j,k,symi,symj,symk,rc1,
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
c
       if (typt3.eq.3) then
c4.4.2add part T2 . U
       call t3sgl(Work(iOff),wrksize,
     & mapdv,symijk,mapdfk4,mapifk4,mapdfk3,mapifk3,
     & mapdt22,mapit22,mapdt21,mapit21,
     & 1,i,j,k,symi,symj,symk,rc1,
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
c
c
c4.5  divide by denominators and calc energy contribution
c
c
c4.5.1divide by den.
       call t3div(Work(iOff),wrksize,
     & mapdw,mapdv,symijk,mapddp2,mapidp2,mapddp1,
     & mapidp1,1,i,j,k,symi,symj,symk,ec,rc1)
c
c4.5.2add energy contribution
       ebbb=ebbb+ec
c
       end if
c
c
 1000   continue
 1100   continue
c
cpar   Separate printing of partial energies e... are
c      useful only in serial run. For paralell run also
c      cycle over k is segmented (via paralelization),
c      so e... are not complete contributions (only sum
c      over all nodes have some sense). Thus, these values
c      in paralell run are too dangerous to use separately
c      so their printout is supressed.
c
       if (nProcs.eq.1) then
         call t3wresult (symi,symj,i,j,eaaa(1),eaab(1),eabb(1),ebbb(1))
         if (fullprint.gt.1) then
           write(6,*) ' Eaaa =',eaaa
           write(6,*) ' Eaab =',eaab
           write(6,*) ' Eabb =',eabb
           write(6,*) ' Ebbb =',ebbb
         end if
       end if
cendpar
c
 1201   continue
 1200   continue
 1300   continue
 1400   continue
c
c
co    ***** final section *****
c
co.*        allreduced energy components
        call gadgop (eaaa,1,'+')
        call gadgop (eaab,1,'+')
        call gadgop (eabb,1,'+')
        call gadgop (ebbb,1,'+')
cstare  call MPI_ALLREDUCE (ebbb,ec,1,
c    c  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)
c       ebbb=ec
c
c
co.*  type results
c
      IF (fullprint.GE.0) THEN
         write(6,'(6X,A,F24.13)')' CCSD     =',eccsd
         write(6,'(6X,A,F24.13)')' T3 corr. =',eaaa+eaab+eabb+ebbb
         write(6,'(6X,A,F24.13)')' CCSD + T3=',eccsd+eaaa+eaab+eabb+ebbb
         write(6,*) ' T3 energy decomposition into spin parts'
         write(6,*) ' Eaaa =',eaaa
         write(6,*) ' Eaab =',eaab
         write(6,*) ' Eabb =',eabb
         write(6,*) ' Ebbb =',ebbb
         write(6,*)
         write(6,*)
         write(6,'(6X,A)') 'Happy Landing!'
         write(6,*)
      ENDIF
      call add_Info('E_CCSD_T',eccsd+eaaa+eaab+eabb+ebbb,1,8)
c Export a method and energy to the MOLCAS runfile
      Call Put_cArray('Relax Method','CCSDT   ',8)
      Call Store_Energies(1,eccsd+eaaa+eaab+eabb+ebbb,1)
c     Releasing the memory
      Call GetMem('CCT3','Free','Real',iOff,wrksize)
c
c
      ireturn=0
      return
      end
c
c     --------------------------------------
c
       subroutine cct3_getint (wrk,wrksize,
     & i,symi,possr0,mapdr,mapir,rc)
c
c     this routine read integrals R_i(a,bc) for given i in given symi
c
c     i      - number of orbital (I)
c     symi   - irrep of i (I)
c     possr0 - initial possition of R (I)
c     mapdr  - direct map of R (I)
c     mapir  - inverse map of R (I)
c     rc     - return (error) code (O)
c
       implicit none
#include "t31.fh"
#include "wrk.fh"
c
       integer i,symi,possr0,rc
c
       integer mapdr(0:512,1:6)
       integer mapir(1:8,1:8,1:8)
c
c     help variables
c
       integer iadd,lun,isym,num
c      integer rc1
       integer poss,lenght,im
c
c1    some tests
c
       if (i.gt.noa(symi)) then
c     RC=1 : i is higher than occipied in this irrep (Stup)
       rc=1
       return
       end if
c
       if (i.lt.1) then
c     RC=2 : i is less than 1 (Stup)
       rc=2
       return
       end if
c
c2    calc number for this orbital
c
       iadd=0
       if (symi.gt.1) then
       do 10 isym=1,symi-1
       iadd=iadd+noa(isym)
 10     continue
       end if
c
       num=iadd+i
c
c
c3    get R
c
       lun=1
       daddr(lun)=T3IntPoss(num)
c
       call daname (lun,t3nam)
c
       call idafile (lun,2,mapdr,513*6,daddr(lun))
       call idafile (lun,2,mapir,8*8*8,daddr(lun))
c
       poss=possr0
       lenght=0
       do im=1,mapdr(0,5)
         mapdr(im,1)=poss
         poss=poss+mapdr(im,2)
         lenght=lenght+mapdr(im,2)
c        write (*,99) ' MAP',(mapdr(im,k),k=1,6)
c99       format (a3,i8,2x,i8,4(2x,i2))
       end do
c
       if (lenght.gt.0) then
         call ddafile (lun,2,wrk(possr0),lenght,daddr(lun))
       end if
c      call cct3_getmediate (wrk,wrksize,
c    & lun,possr0,mapdr,mapir,rc1)

       call daclos (lun)
c
       return
       end
c
c     --------------------------------------
c
       subroutine GetIntPoss
c
c        this routine read T3IntPoss array from the first record
c        of t3nam file
c
        implicit none
#include "t31.fh"
c
        integer lun
c
       lun=1
       call daname (lun,t3nam)
       daddr(lun)=0
       call idafile (lun,2,T3IntPoss,maxorb,daddr(lun))
       call daclos (lun)
c
       return
       end
c
c     ----------------------------------------
c
