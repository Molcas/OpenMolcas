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
       subroutine CCT3(ireturn)
!
!     this program calculate noniterative T3 contributions
!     to CCSD energy
!
      use Para_Info, only: MyRank, nProcs
#include "t31.fh"
#include "t32.fh"
!
!     work file declaration
       integer wrksize
#include "WrkSpc.fh"
!
!     variables
!
       real*8 eccsd,eaaa(1),eaab(1),eabb(1),ebbb(1),ec
       integer i,j,k,i3,j3,k3,symi,symj,symk,symi3,symj3,symk3
       integer symijk,symij
       integer nsg,keyyes,rc1,jup,ilow,posst,ssl1,ssl2,ssm1,ssm2
       integer symjstart,symjstop
       integer jstart,jstop,istart,istop
!par
       integer id,counter
!
!       integer nhelp
!
!
!par
!stare Call SetTim
!stare call MPI_COMM_RANK(MPI_COMM_WORLD,myRank,rc)
!stare call MPI_COMM_SIZE(MPI_COMM_WORLD,nProcs,rc)

!I    **************  start   section **************
!
!I.1  read input data form INPDAT (reorg) and input file
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
!
!I.2  write head to output file
      If (fullprint.GE.0) call t3wrhead
!
!I.3.1 calc. work space requirements to fix and help mediates
      call t3initfiles (wrksize)
      If (fullprint.GE.0)                                               &
     & write(6,*) ' Work space requirements :',wrksize
!
!I.3.2 allocate work space
!
      Call GetMem('CCT3','Max','Real',maxspace,maxspace)
      maxspace=maxspace-4
      if (maxspace.lt.wrksize) then
         write(6,*) ' Allocation of work space failed',maxspace
         write(6,*) ' Increase the value of the variable MOLCAS_MEM'
         Call Abend()
      end if
      Call GetMem('CCT3','Allo','Real',iOff,wrksize)
!
!I.3.3 set wrk = 0
      call cct3_mv0zero (wrksize,wrksize,Work(iOff))
      If (fullprint.GE.0) write(6,*) ' Allocation of work space : Done'
!
!I.4  read static integrals from INTSTA (reorg) file
      call t3reaintsta (Work(iOff),wrksize)
!
!I.5  divide fok to faa,fai,fii and dp
      call cct3_divfok (Work(iOff),wrksize,                             &
     & mapdn,mapin,possn0,mapdp,mapip,possp0,                           &
     & mapdfk1,mapifk1,possfk10,mapdfk2,mapifk2,possfk20,               &
     & mapdfk3,mapifk3,possfk30,mapdfk4,mapifk4,possfk40,               &
     & mapdfk5,mapifk5,possfk50,mapdfk6,mapifk6,possfk60,               &
     & mapddp1,mapidp1,possdp10,mapddp2,mapidp2,possdp20,rc1)
!
!I.6  read CCSD energy and amplitudes
      call t3reaccsd (Work(iOff),wrksize,                               &
     & eccsd)
!
!I.7  get address vector T3IntPos
!       they are located at the beggining of the t3nam file
      call GetIntPoss
!
!I.*  set energies=0
      eaaa=0.0d0
      eaab=0.0d0
      eabb=0.0d0
      ebbb=0.0d0
!
!I.par initialize parallel counter
      counter=0
!
!
!     ***** work section *****
!
!
      if (ijsegkey.eq.0) then
!No Segmentation
         symimin=1
         symimax=nsym
      end if
!
!Segmented
!noseg do 1400 symi=1,nsym
      do 1400 symi=symimin,symimax
       if (noa(symi).eq.0) then
       goto 1400
       end if
!
!      def symjstart,symjstop
!
       if (ijsegkey.eq.0) then
!No Segmentation
          symjstart=1
          symjstop=symi
       else
!
!Segmented
       if (symimin.eq.symjmax) then
!      case symimin=symi=symimax
         symjstart=symjmin
         symjstop=symjmax
       else
!      case symimin<symjmax
!
         if(symi.eq.symimin) then
           symjstart=symjmin
           symjstop=symi
         else if (symi.eq.symimax) then
           symjstart=1
           symjstop=symjmax
         else
!c       sub case symimin<symi<symimax
           symjstart=1
           symjstop=symi
         end if
!
       end if
       end if
!
!
!noseg do 1300 symj=1,symi
       do 1300 symj=symjstart,symjstop
       if ((symi.eq.symj).and.(noa(symi).le.1)) then
       goto 1300
       else if (noa(symj).eq.0) then
       goto 1300
       end if
!
!*    define sym(ij)
       symij=mmul(symi,symj)
!
!*    write symmetry starus of I and J
       if (fullprint.gt.0) then
       write(6,*) ' SYMI',symi,'SYMJ',symj
       end if
!
!
       if (symi.eq.symj) then
       ilow=2
       else
       ilow=1
       end if
!
!
!*    loop over inex I
!
!       define istart,istop
!
       if (ijsegkey.eq.0) then
!No Segmentation
          istart=ilow
          istop=noa(symi)
       else
!
!Segmented
        if ((symi.eq.symimin).and.(symj.eq.symjmin)) then
          istart=imin
        else
          istart=ilow
        end if
!
        if ((symi.eq.symimax).and.(symj.eq.symjmax)) then
          istop=imax
        else
          istop=noa(symi)
        end if
!
       end if
!
!noseg    do 1200 i=ilow,noa(symi)
          do 1200 i=istart,istop
!
!*    get integrals <ab|ic> for given i into R1(a,bc)
       call cct3_getint (Work(iOff),wrksize,                            &
     & i,symi,possr10,mapdr1,mapir1,rc1)
!
!*    def upper limit for index j
       if (symi.eq.symj) then
       jup=i-1
       else
       jup=noa(symj)
       end if
!
!
!*    loop over index J
!
!*      def jstart,jstop
!
       if (ijsegkey.eq.0) then
!No Segmentation
          jstart=1
          jstop=jup
       else
!
!Segmented
         if ((symimin.eq.symimax).and.(symjmin.eq.symjmax)) then
!        we are in the only symmetry block taken
!        into the consideration
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
!        we are in initial symmetry block
           if (i.eq.imin) then
             jstart=jmin
             jstop=jup
           else
             jstart=1
             jstop=jup
           end if
         else if ((symi.eq.symimax).and.(symj.eq.symjmax)) then
!        we are in terminal symmetry block
           if (i.eq.imax) then
             jstart=1
             jstop=jmax
           else
             jstart=1
             jstop=jup
           end if
         else
!        we are in intermediate symmetry block
           jstart=1
           jstop=jup
         end if
!
       end if
!
!
!noseg do 1200 j=1,jup
       do 1201 j=jstart,jstop
!
!*    get integrals <ab|jc> for given j into R2(a,bc)
       call cct3_getint (Work(iOff),wrksize,                            &
     & j,symj,possr20,mapdr2,mapir2,rc1)
!
!
       do 1100 symk=1,nsym
       if (noa(symk).eq.0) then
       goto 1100
       end if
!
       if (fullprint.gt.1) then
       write(6,*) ' SYMI',symi,'SYMJ',symj,'SYMK',symk
       end if
!
!*    define sym(ijk)
       symijk=mmul(symij,symk)
!
!
!*    loop over index K
!
       do 1000 k=1,noa(symk)
       if (fullprint.ge.2) then
       write (6,999) i,j,k
999    format (' I,J,K ',3(i3,1x))
       end if
!
!par    update parallel counter, choose proper id for this 'portion'
!       and skip if this portion is not for myRank
        counter=counter+1
        id=mod(counter,nProcs)
        if (myRank.ne.id) goto 1000
!
!*    get integrals <ab|kc> for given k into R3(a,bc)
       call cct3_getint (Work(iOff),wrksize,                            &
     & k,symk,possr30,mapdr3,mapir3,rc1)
!
!
!1    ***** aaa spin combination *****
!
!
!1*   def keyyes
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
!
       if (keyyes.eq.1) then
!
!1.*  define maps of W(abc)
       call cct3_t3grc0(3,5,3,3,3,0,symijk,possw0,posst,mapdw,mapiw)
!
!1.*  vanish W(abc)
       call stz (Work(iOff),wrksize,                                    &
     & mapdw)
!
!
!1.1  permutations (ijk) P(a,bc) (general sign+)
       nsg=1
!
!
!1.1.1V graph
!
!1.1.1*            def L1(bc,d) <- R3(b,cd) for given k
       call defv (Work(iOff),wrksize,                                   &
     & 1,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr3,mapir3,symk,rc1)
!
!1.1.1*            ext M1(da) <- T2aaaa(da,ij) for given i,j
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,j,0,symi,symj,0,mapdt21,mapit21,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!1.1.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,                            &
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,                    &
     & rc1)
       ssm2=ssm1
!
!1.1.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!1.1.1*            pack W(abc) <-  P(a,bc) [L2(bc,a)] (minus due to the
!     usung Tikda instead of T2ikad)
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
!
!1.1.2O graph
!
!1.1.2*            ext L1(bc,l) <- T2aaaa(bc,kl) for given k
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,k,0,0,symk,0,0,mapdt21,mapit21,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!1.1.2*     ext M1(l,a) <- W11(l,a,ij)=<la||ij>aaaa for given ij
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,j,0,symi,symj,0,mapdw11,mapiw11,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!1.1.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!1.1.2*            pack W(abc) <- P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!
!1.2  permutations (ikj) P(a,bc)
       nsg=-1
!
!
!1.2.1V graph
!
!1.2.1*            def L1(bc,d) <- R2(b,cd) for given j
       call defv (Work(iOff),wrksize,                                   &
     & 1,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr2,mapir2,symj,rc1)
!
!1.2.1*            ext M1(da) <- T2aaaa(da,ik) for given i,k
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,k,0,symi,symk,0,mapdt21,mapit21,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!1.2.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,                            &
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,                    &
     & rc1)
       ssm2=ssm1
!
!1.2.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!1.2.1*            pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to
!     usung Tikda instead of T2ikad)
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
!
!1.2.2O graph
!
!1.2.2*            ext L1(bc,l) <- T2aaaa(bc,jl) for given j
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,j,0,0,symj,0,0,mapdt21,mapit21,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!1.2.2*     ext M1(l,a) <- W11(l,a,ik)=<la||ik>aaaa for given ik
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,k,0,symi,symk,0,mapdw11,mapiw11,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!1.2.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!1.2.2*            pack W(abc) <- P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!
!1.3  permutations (jki) P(a,bc)
       nsg=1
!
!
!1.3.1V graph
!
!1.3.1*            def L1(bc,d) <- R1(b,cd) for given i
       call defv (Work(iOff),wrksize,                                   &
     & 1,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr1,mapir1,symi,rc1)
!
!1.3.1*            ext M1(da) <- T2aaaa(da,jk) for given j,k
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j,k,0,symj,symk,0,mapdt21,mapit21,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!1.3.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,                            &
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,                    &
     & rc1)
       ssm2=ssm1
!
!1.3.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!1.3.1*            pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to
!     usung Tjkda instead of T2jkad)
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
!
!1.3.2O graph
!
!1.3.2*            ext L1(bc,l) <- T2aaaa(bc,kl) for given i
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,i,0,0,symi,0,0,mapdt21,mapit21,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!1.3.2*     ext M1(l,a) <- W11(l,a,jk)=<la||jk>aaaa for given jk
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j,k,0,symj,symk,0,mapdw11,mapiw11,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!1.3.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!1.3.2*            pack W(abc) <-  P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!
!1.4  add singles
!
!1.4.0mov V <- -W
       call cct3_t3grc0(3,5,3,3,3,0,symijk,possv0,posst,mapdv,mapiv)
       call minusa (Work(iOff),wrksize,                                 &
     & mapdw,-1.0d0)
       call setb (Work(iOff),wrksize,                                   &
     & mapdw,mapdv,1.0d0)
!
       if (typt3.gt.1) then
!1.4.1add part W2 . T1
       call t3sgl(Work(iOff),wrksize,                                   &
     & mapdv,symijk,mapdt11,mapit11,mapdt12,mapit12,                    &
     & mapdw21,mapiw21,mapdw22,mapiw22,                                 &
     & 1,i,j,k,symi,symj,symk,rc1,                                      &
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,                     &
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,                     &
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
!
       if (typt3.eq.3) then
!1.4.2add part T2 . U
       call t3sgl(Work(iOff),wrksize,                                   &
     & mapdv,symijk,mapdfk3,mapifk3,mapdfk4,mapifk4,                    &
     & mapdt21,mapit21,mapdt22,mapit22,                                 &
     & 1,i,j,k,symi,symj,symk,rc1,                                      &
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,                     &
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,                     &
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
!
!
!1.5  divide by denominators and calc energy contribution
!
!
!1.5.1divide by den.
       call t3div(Work(iOff),wrksize,                                   &
     & mapdw,mapdv,symijk,mapddp1,mapidp1,mapddp2,                      &
     & mapidp2,1,i,j,k,symi,symj,symk,ec,rc1)
!
!1.5.2add energy contribution
       eaaa=eaaa+ec
!
       end if
!
!
!2    ***** aab spin combination *****
!
!
!2*   def keyyes
       if (k.le.nob(symk)) then
       keyyes=1
       else
       keyyes=0
       end if
!
       if (keyyes.eq.1) then
!
!2.*  define maps of W(abc)
       call cct3_t3grc0(3,1,3,3,4,0,symijk,possw0,posst,mapdw,mapiw)
!
!2.*  vanish W(abc)
       call stz (Work(iOff),wrksize,                                    &
     & mapdw)
!
!
!2.1  permutations (ijk) P(a,b) (c)
       nsg=1
!
!2.1.1V graph
!
!2.1.1*            def L1(b,c,d)aba <- R3(b,cd) for given k
       call defv (Work(iOff),wrksize,                                   &
     & 4,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr3,mapir3,symk,rc1)
!
!2.1.1*            ext M1(da) <- T2aaaa(da,ij) for given i,j
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,j,0,symi,symj,0,mapdt21,mapit21,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.1.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,                            &
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,                    &
     & rc1)
       ssm2=ssm1
!
!2.1.1*            mult L2(b,c,a) <- L1(b,c,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.1.1*            pack W(ab,c) <- - P(a,b) [L2(b,c,a)] (minus is due to
!     usung Tijda instead of T2ijad)
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
!
!2.1.2O graph
!
!2.1.2*            ext L1(b,c,l) <- T2abab(b,c,l,k) for given k
       call ext(Work(iOff),wrksize,                                     &
     & 4,4,k,0,0,symk,0,0,mapdt23,mapit23,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!2.1.2*     ext M1(l,a) <- W12(l,a,ij)=<la||ij>aaaa for given ij
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,j,0,symi,symj,0,mapdw11,mapiw11,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.1.2*            mult L2(b,c,a) <- L1(b,c,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.1.2*            pack W(ab,c) <- - P(a,b) [L2(b,c,a)] (-, premutation in V)
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
!
!
!2.2  permutations (ijk) (cab) do not nontribute
!
!
!2.3  permutations (ikj) P(a,b) (c)
       nsg=-1
!
!2.3.1V graph
!
!2.3.1*            def L1(b,c,d)abb <- R2(b,cd) for given j
       call defv (Work(iOff),wrksize,                                   &
     & 3,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr2,mapir2,symj,rc1)
!
!2.3.1*            ext M1(a,d) <- T2abab(a,d,i,k) for given i,k
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,k,0,symi,symk,0,mapdt23,mapit23,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.3.1*            map M2(d,a) <- M1(a,d)
       call cct3_map (Work(iOff),wrksize,                               &
     & 2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,                      &
     & possm20,posst,rc1)
       ssm2=ssm1
!
!2.3.1*            mult L2(b,c,a) <- L1(b,c,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.3.1*            pack W(ab,c) <-  P(a,b) [L2(b,c,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!2.3.2O graph
!
!2.3.2*            ext L1(b,c,l) <- T2abab(b,c,j,l) for given j
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,j,0,0,symj,0,0,mapdt23,mapit23,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!2.3.2*     ext M1(l,a) <- W14(l,a,ik)=<la||ik>baab for given ik
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,k,0,symi,symk,0,mapdw14,mapiw14,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.3.2*            mult L2(b,c,a) <- L1(b,c,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.3.2*            pack W(ab,c) <- P(a,b) [L2(b,c,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!
!2.4  permutations (ikj) (cab)
       nsg=-1
!
!
!2.4.1V graph
!
!2.4.1*            def L1(ab,d)aaa <- R2(a,bc) for given j
       call defv (Work(iOff),wrksize,                                   &
     & 1,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr2,mapir2,symj,rc1)
!
!2.4.1*            ext M1(d,c) <- T2abab(d,c,i,k) for given i,k
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,k,0,symi,symk,0,mapdt23,mapit23,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.4.1*            mult L2(ab,c) <- L1(ab,d) . M1(d,c)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.4.1*            add (pack) W(ab,c) <-  [L2(b,c,a)] (- due to permuted T)
       call cct3_add (Work(iOff),wrksize,                               &
     & 3,3,0,0,0,0,1,1,-1.0d0*nsg,mapdl2,symijk,                        &
     & mapdw,mapiw,symijk,rc1)
!
!2.4.2O graph
!
!2.4.2*            ext L1(ab,l) <- T2aaaa(ab,jl) for given j
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,j,0,0,symj,0,0,mapdt21,mapit21,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!2.4.2*     ext M1(l,c) <- W13(l,c,ik)=<lc||ik>abab for given ik
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,k,0,symi,symk,0,mapdw13,mapiw13,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.4.2*            mult L2(ab,c) <- L1(ab,l) . M1(l,c)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.4.2*            add (pack) W(ab,c) <- [L2(ab,c)]
       call cct3_add (Work(iOff),wrksize,                               &
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl2,symijk,                         &
     & mapdw,mapiw,symijk,rc1)
!
!
!2.5  permutations (jki) P(a,b) (c)
       nsg=1
!
!
!2.5.1V graph
!
!2.5.1*            def L1(b,c,d)abb <- R1(b,cd) for given i
       call defv (Work(iOff),wrksize,                                   &
     & 3,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr1,mapir1,symi,rc1)
!
!2.5.1*            ext M1(a,d) <- T2abab(a,d,j,k) for given j,k
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j,k,0,symj,symk,0,mapdt23,mapit23,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.5.1*            map M2(d,a) <- M1(a,d)
       call cct3_map (Work(iOff),wrksize,                               &
     & 2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,                      &
     & possm20,posst,rc1)
       ssm2=ssm1
!
!2.5.1*            mult L2(b,c,a) <- L1(b,c,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.5.1*            pack W(ab,c) <-  P(a,b) [L2(b,c,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!2.5.2O graph
!
!2.5.2*            ext L1(b,c,l) <- T2abab(b,c,i,l) for given i
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,i,0,0,symi,0,0,mapdt23,mapit23,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!2.5.2*     ext M1(l,a) <- W14(l,a,jk)=<la||jk>baab for given jk
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j,k,0,symj,symk,0,mapdw14,mapiw14,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.5.2*            mult L2(b,c,a) <- L1(b,c,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.5.2*            pack W(ab,c) <- P(a,b) [L2(b,c,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,2,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!
!2.6  permutations (ikj) (cab)
       nsg=1
!
!
!2.6.1V graph
!
!2.6.1*            def L1(ab,d)aaa <- R1(a,bc) for given i
       call defv (Work(iOff),wrksize,                                   &
     & 1,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr1,mapir1,symi,rc1)
!
!2.6.1*            ext M1(d,c) <- T2abab(d,c,j,k) for given j,k
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j,k,0,symj,symk,0,mapdt23,mapit23,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.6.1*            mult L2(ab,c) <- L1(ab,d) . M1(d,c)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.6.1*            add (pack) W(ab,c) <-  [L2(b,c,a)] (- due to permuted T)
       call cct3_add (Work(iOff),wrksize,                               &
     & 3,3,0,0,0,0,1,1,-1.0d0*nsg,mapdl2,symijk,                        &
     & mapdw,mapiw,symijk,rc1)
!
!2.6.2O graph
!
!2.6.2*            ext L1(ab,l) <- T2aaaa(ab,il) for given i
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,i,0,0,symi,0,0,mapdt21,mapit21,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!2.6.2*     ext M1(l,c) <- W13(l,c,jk)=<lc||jk>abab for given jk
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j,k,0,symj,symk,0,mapdw13,mapiw13,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!2.6.2*            mult L2(ab,c) <- L1(ab,l) . M1(l,c)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!2.6.2*            add (pack) W(ab,c) <- [L2(ab,c)]
       call cct3_add (Work(iOff),wrksize,                               &
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl2,symijk,                         &
     & mapdw,mapiw,symijk,rc1)
!
!
!2.7  add singles
!
!
!2.7.0mov V <- W
       call cct3_t3grc0(3,1,3,3,4,0,symijk,possv0,posst,mapdv,mapiv)
       call minusa (Work(iOff),wrksize,                                 &
     & mapdw,-1.0d0)
       call setb (Work(iOff),wrksize,                                   &
     & mapdw,mapdv,1.0d0)
!
       if (typt3.gt.1) then
!2.7.1add part W2 . T1
       call t3sgl(Work(iOff),wrksize,                                   &
     & mapdv,symijk,mapdt11,mapit11,mapdt12,mapit12,                    &
     & mapdw21,mapiw21,mapdw23,mapiw23,                                 &
     & 2,i,j,k,symi,symj,symk,rc1,                                      &
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,                     &
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,                     &
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
!
       if (typt3.eq.3) then
!2.7.2add part T2 . U
       call t3sgl(Work(iOff),wrksize,                                   &
     & mapdv,symijk,mapdfk3,mapifk3,mapdfk4,mapifk4,                    &
     & mapdt21,mapit21,mapdt23,mapit23,                                 &
     & 2,i,j,k,symi,symj,symk,rc1,                                      &
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,                     &
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,                     &
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
!
!
!2.8  divide by denominators and calc energy contribution
!
!
!2.8.1divide by den.
       call t3div(Work(iOff),wrksize,                                   &
     & mapdw,mapdv,symijk,mapddp1,mapidp1,mapddp2,                      &
     & mapidp2,2,i,j,k,symi,symj,symk,ec,rc1)
!
!2.8.2add energy contribution
       eaab=eaab+ec
!
       end if
!
!
!3    ***** abb spin combination *****
!
!
!     Note:
!     in spin combination 3-abb indexes are changed as follows
!
       i3=k
       symi3=symk
       j3=i
       symj3=symi
       k3=j
       symk3=symj
!
!     therefore, also R files are mixed:
!     R13 = R3, R23=R1 and R33=R2, symijk3=symijk
!
!3.*  def keyyes
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
!
       if (keyyes.eq.1) then
!
!3.*  define maps of W(a,bc)
       call cct3_t3grc0(3,2,3,4,4,0,symijk,possw0,posst,mapdw,mapiw)
!
!3.*  vanish W(abc)
       call stz (Work(iOff),wrksize,                                    &
     & mapdw)
!
!
!3.1  permutation (ijk) (abc)
       nsg=1
!
!
!3.1.1V graph
!
!3.1.1*            def L1(bc,d) <- R33(b,cd) for given k3
       call defv (Work(iOff),wrksize,                                   &
     & 2,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr2,mapir2,symk3,rc1)
!
!3.1.1*            ext M1(a,d) <- T2abab(a,d,i3,j3) for given i3,j3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i3,j3,0,symi3,symj3,0,mapdt23,mapit23,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.1.1*            map M2(d,a) <- M1(a,d)
       call cct3_map (Work(iOff),wrksize,                               &
     & 2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,                      &
     & possm20,posst,rc1)
       ssm2=ssm1
!
!3.1.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.1.1map L1(a,bc) <- L2(bc,a)
       call cct3_map (Work(iOff),wrksize,                               &
     & 3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,                      &
     & possl10,posst,rc1)
       ssl1=ssl2
!3.1.1*            add (pack) W(a,bc) <- [L1(a,bc)]
       call cct3_add (Work(iOff),wrksize,                               &
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl1,symijk,                         &
     & mapdw,mapiw,symijk,rc1)
!
!3.1.2O graph
!
!3.1.2*            ext L1(bc,l) <- T2bbbb(bc,k3l) for given k3
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,k3,0,0,symk3,0,0,mapdt22,mapit22,1,                          &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!3.1.2*     ext M1(l,a) <- W14(l,a,i3j3)=<la||i3j3>baab for given i3,j3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i3,j3,0,symi3,symj3,0,mapdw14,mapiw14,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.1.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.1.2map L1(a,bc) <- L2(bc,a)
       call cct3_map (Work(iOff),wrksize,                               &
     & 3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,                      &
     & possl10,posst,rc1)
       ssl1=ssl2

!3.1.2*            add (pack) W(a,bc) <- [L1(a,bc)]
       call cct3_add (Work(iOff),wrksize,                               &
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl1,symijk,                         &
     & mapdw,mapiw,symijk,rc1)
!
!
!3.2  permutations (ijk)(bac),(ijk)(cab)
       nsg=1
!
!
!3.2.1V graph
!
!3.2.1*            def L1(a,c,d)aba <- R33(a,cd) for given k3
       call defv (Work(iOff),wrksize,                                   &
     & 4,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr2,mapir2,symk3,rc1)
!
!3.2.1*            ext M1(d,b) <- T2abab(d,b,i3,j3) for given i3,j3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i3,j3,0,symi3,symj3,0,mapdt23,mapit23,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.2.1*            mult L2(ac,b) <- L1(a,c,d) . M1(d,b)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.2.1*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm T)
       call t3addpck (Work(iOff),wrksize,                               &
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,                         &
     & 0,rc1)
!
!3.2.2O graph
!
!3.2.2*            ext L1(a,c,l) <- T2abab(a,c,l,k3) for given k3
       call ext(Work(iOff),wrksize,                                     &
     & 4,4,k3,0,0,symk3,0,0,mapdt23,mapit23,1,                          &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!3.2.2*     ext M1(l,b) <- W13(l,b,i3j3)=<la||i3j3>abab for given i3,j3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i3,j3,0,symi3,symj3,0,mapdw13,mapiw13,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.2.2*            mult L2(ac,b) <- L1(a,c,l) . M1(l,b)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.2.2*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm V)
       call t3addpck (Work(iOff),wrksize,                               &
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,                         &
     & 0,rc1)
!
!
!3.3  permutation (ikj) (abc)
       nsg=-1
!
!
!3.3.1V graph
!
!3.3.1*            def L1(bc,d) <- R23(b,cd) for given j3
       call defv (Work(iOff),wrksize,                                   &
     & 2,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr1,mapir1,symj3,rc1)
!
!3.3.1*            ext M1(a,d) <- T2abab(a,d,i3,k3) for given i3,k3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i3,k3,0,symi3,symk3,0,mapdt23,mapit23,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.3.1*            map M2(d,a) <- M1(a,d)
       call cct3_map (Work(iOff),wrksize,                               &
     & 2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,                      &
     & possm20,posst,rc1)
       ssm2=ssm1
!
!3.3.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.3.1map L1(a,bc) <- L2(bc,a)
       call cct3_map (Work(iOff),wrksize,                               &
     & 3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,                      &
     & possl10,posst,rc1)
       ssl1=ssl2

!3.3.1*            add (pack) W(a,bc) <- [L1(a,bc)]
       call cct3_add (Work(iOff),wrksize,                               &
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl1,symijk,                         &
     & mapdw,mapiw,symijk,rc1)
!
!3.3.2O graph
!
!3.3.2*            ext L1(bc,l) <- T2bbbb(bc,j3l) for given j3
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,j3,0,0,symj3,0,0,mapdt22,mapit22,1,                          &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!3.3.2*     ext M1(l,a) <- W14(l,a,i3k3)=<la||i3k3>baab for given i3,k3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i3,k3,0,symi3,symk3,0,mapdw14,mapiw14,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.3.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.3.2map L1(a,bc) <- L2(bc,a)
       call cct3_map (Work(iOff),wrksize,                               &
     & 3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,                      &
     & possl10,posst,rc1)
       ssl1=ssl2

!3.3.2*            add (pack) W(a,bc) <- [L1(a,bc)]
       call cct3_add (Work(iOff),wrksize,                               &
     & 3,3,0,0,0,0,1,1,1.0d0*nsg,mapdl1,symijk,                         &
     & mapdw,mapiw,symijk,rc1)
!
!
!3.4  permutations (ikj)(bac),(ikj)(cab)
       nsg=-1
!
!
!3.4.1V graph
!
!3.4.1*            def L1(a,c,d)aba <- R23(a,cd) for given j3
       call defv (Work(iOff),wrksize,                                   &
     & 4,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr1,mapir1,symj3,rc1)
!
!3.4.1*            ext M1(d,b) <- T2abab(d,b,i3,k3) for given i3,k3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i3,k3,0,symi3,symk3,0,mapdt23,mapit23,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.4.1*            mult L2(ac,b) <- L1(a,c,d) . M1(d,b)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.4.1*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm T)
       call t3addpck (Work(iOff),wrksize,                               &
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,                         &
     & 0,rc1)
!
!3.4.2O graph
!
!3.4.2*            ext L1(a,c,l) <- T2abab(a,c,l,j3) for given j3
       call ext(Work(iOff),wrksize,                                     &
     & 4,4,j3,0,0,symj3,0,0,mapdt23,mapit23,1,                          &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!3.4.2*     ext M1(l,b) <- W13(l,b,i3k3)=<la||i3k3>abab for given i3,k3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i3,k3,0,symi3,symk3,0,mapdw13,mapiw13,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.4.2*            mult L2(ac,b) <- L1(a,c,l) . M1(l,b)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.4.2*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm V)
       call t3addpck (Work(iOff),wrksize,                               &
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,                         &
     & 0,rc1)
!
!
!3.5  permutations (jki)(abc) do not contribute
!
!
!3.6  permutations (jki)(bac),(jki)(cab)
       nsg=1
!
!
!3.6.1V graph
!
!3.6.1*            def L1(a,c,d)abb <- R13(a,cd) for given i3
       call defv (Work(iOff),wrksize,                                   &
     & 3,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr3,mapir3,symi3,rc1)
!
!3.6.1*            ext M1(db) <- T2bbbb(db,j3k3) for given j3,k3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j3,k3,0,symj3,symk3,0,mapdt22,mapit22,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.6.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,                            &
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,                    &
     & rc1)
       ssm2=ssm1
!
!3.6.1*            mult L2(ac,b) <- L1(a,c,d) . M2(d,b)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.6.1*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm T)
       call t3addpck (Work(iOff),wrksize,                               &
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,                         &
     & 0,rc1)
!
!3.6.2O graph
!
!3.6.2*            ext L1(a,c,l) <- T2abab(a,c,i3,l) for given i3
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,i3,0,0,symi3,0,0,mapdt23,mapit23,1,                          &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!3.6.2*     ext M1(l,b) <- W12(l,b,j3k3)=<la||j3k3>bbbb for given j3,k3
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j3,k3,0,symj3,symk3,0,mapdw12,mapiw12,1,                     &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!3.6.2*            mult L2(ac,b) <- L1(a,c,l) . M1(l,b)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!3.6.2*            add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)]
       call t3addpck (Work(iOff),wrksize,                               &
     & 3,3,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,                          &
     & 0,rc1)
!
!
!3.7  add singles
!
!
!3.7.0mov V <- W
       call cct3_t3grc0(3,2,3,4,4,0,symijk,possv0,posst,mapdv,mapiv)
       call minusa (Work(iOff),wrksize,                                 &
     & mapdw,-1.0d0)
       call setb (Work(iOff),wrksize,                                   &
     & mapdw,mapdv,1.0d0)
!
       if (typt3.gt.1) then
!3.7.1add part W2 . T1
       call t3sgl(Work(iOff),wrksize,                                   &
     & mapdv,symijk,mapdt11,mapit11,mapdt12,mapit12,                    &
     & mapdw23,mapiw23,mapdw22,mapiw22,                                 &
     & 3,i3,j3,k3,symi3,symj3,symk3,rc1,                                &
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,                     &
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,                     &
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
!
       if (typt3.eq.3) then
!3.7.2add part T2 . U
       call t3sgl(Work(iOff),wrksize,                                   &
     & mapdv,symijk,mapdfk3,mapifk3,mapdfk4,mapifk4,                    &
     & mapdt23,mapit23,mapdt22,mapit22,                                 &
     & 3,i3,j3,k3,symi3,symj3,symk3,rc1,                                &
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,                     &
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,                     &
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
!
!
!3.8  divide by denominators and calc energy contribution
!
!
!3.8.1divide by den.
       call t3div(Work(iOff),wrksize,                                   &
     & mapdw,mapdv,symijk,mapddp1,mapidp1,mapddp2,                      &
     & mapidp2,3,i3,j3,k3,symi3,symj3,symk3,ec,rc1)
!
!3.8.2add energy contribution
       eabb=eabb+ec
!
       end if
!
!
!4    ***** bbb spin combination *****
!
!
!4.*  def keyyes
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
!
       if ((i.gt.nob(symi)).or.(j.gt.nob(symj))                         &
     & .or.(k.gt.nob(symk))) then
       keyyes=0
       end if
!
       if (keyyes.eq.1) then
!
!4.*  define maps of W(abc)
       call cct3_t3grc0(3,5,4,4,4,0,symijk,possw0,posst,mapdw,mapiw)
!
!4.*  vanish W(abc)
       call stz (Work(iOff),wrksize,                                    &
     & mapdw)
!
!
!4.1  permutations (ijk) P(a,bc) (general sign+)
       nsg=1
!
!
!4.1.1V graph
!
!4.1.1*            def L1(bc,d) <- R3(b,cd) for given k
       call defv (Work(iOff),wrksize,                                   &
     & 2,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr3,mapir3,symk,rc1)
!
!4.1.1*            ext M1(da) <- T2bbbb(da,ij) for given i,j
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,j,0,symi,symj,0,mapdt22,mapit22,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!4.1.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,                            &
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,                    &
     & rc1)
       ssm2=ssm1
!
!4.1.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!4.1.1*            pack W(abc) <-  P(a,bc) [L2(bc,a)] (minus due to the
!     usung Tikda instead of T2ikad)
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
!
!4.1.2O graph
!
!4.1.2*            ext L1(bc,l) <- T2bbbb(bc,kl) for given k
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,k,0,0,symk,0,0,mapdt22,mapit22,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!4.1.2*     ext M1(l,a) <- W12(l,a,ij)=<la||ij>bbbb for given ij
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,j,0,symi,symj,0,mapdw12,mapiw12,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!4.1.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!4.1.2*            pack W(abc) <- P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!
!4.2  permutations (ikj) P(a,bc)
       nsg=-1
!
!
!4.2.1V graph
!
!4.2.1*            def L1(bc,d) <- R2(b,cd) for given j
       call defv (Work(iOff),wrksize,                                   &
     & 2,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr2,mapir2,symj,rc1)
!
!4.2.1*            ext M1(da) <- T2bbbb(da,ik) for given i,k
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,k,0,symi,symk,0,mapdt22,mapit22,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!4.2.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,                            &
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,                    &
     & rc1)
       ssm2=ssm1
!
!4.2.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!4.2.1*            pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to
!     usung Tikda instead of T2ikad)
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
!
!4.2.2O graph
!
!4.2.2*            ext L1(bc,l) <- T2bbbb(bc,jl) for given j
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,j,0,0,symj,0,0,mapdt22,mapit22,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!4.2.2*     ext M1(l,a) <- W12(l,a,ik)=<la||ik>bbbb for given ik
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,i,k,0,symi,symk,0,mapdw12,mapiw12,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!4.2.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!4.2.2*            pack W(abc) <- P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!
!4.3  permutations (jki) P(a,bc)
       nsg=1
!
!
!4.3.1V graph
!
!4.3.1*            def L1(bc,d) <- R1(b,cd) for given i
       call defv (Work(iOff),wrksize,                                   &
     & 2,possl10,mapdl1,mapil1,ssl1,                                    &
     & mapdr1,mapir1,symi,rc1)
!
!4.3.1*            ext M1(da) <- T2bbbb(da,jk) for given j,k
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j,k,0,symj,symk,0,mapdt22,mapit22,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!4.3.1*            exp M2(d,a) <- M1(da)
       call cct3_expand (Work(iOff),wrksize,                            &
     & 2,1,mapdm1,mapim1,ssm1,possm20,mapdm2,mapim2,                    &
     & rc1)
       ssm2=ssm1
!
!4.3.1*            mult L2(bc,a) <- L1(bc,d) . M2(d,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!4.3.1*            pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to
!     usung Tjkda instead of T2jkad)
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,-nsg,0,rc1)
!
!4.3.2O graph
!
!4.3.2*            ext L1(bc,l) <- T2bbbb(bc,kl) for given i
       call ext(Work(iOff),wrksize,                                     &
     & 4,3,i,0,0,symi,0,0,mapdt22,mapit22,1,                            &
     & possl10,mapdl1,mapil1,ssl1,rc1)
!
!4.3.2*     ext M1(l,a) <- W12(l,a,jk)=<la||jk>bbbb for given jk
       call ext(Work(iOff),wrksize,                                     &
     & 4,7,j,k,0,symj,symk,0,mapdw12,mapiw12,1,                         &
     & possm10,mapdm1,mapim1,ssm1,rc1)
!
!4.3.2*            mult L2(bc,a) <- L1(bc,l) . M1(l,a)
       call cct3_mult(Work(iOff),wrksize,                               &
     & 3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,                   &
     & mapdl2,mapil2,ssl2,possl20,rc1)
!
!4.3.2*            pack W(abc) <-  P(a,bc) [L2(bc,a)]
       call t3addpck(Work(iOff),wrksize,                                &
     & 3,1,mapdl2,mapil2,ssl2,mapdw,mapiw,nsg,0,rc1)
!
!
!4.4  add singles
!
!
!4.4.0mov V <- W
       call cct3_t3grc0(3,5,4,4,4,0,symijk,possv0,posst,mapdv,mapiv)
       call minusa (Work(iOff),wrksize,                                 &
     & mapdw,-1.0d0)
       call setb (Work(iOff),wrksize,                                   &
     & mapdw,mapdv,1.0d0)
!
       if (typt3.gt.1) then
!4.4.1add part W2 . T1
       call t3sgl(Work(iOff),wrksize,                                   &
     & mapdv,symijk,mapdt12,mapit12,mapdt11,mapit11,                    &
     & mapdw22,mapiw22,mapdw21,mapiw21,                                 &
     & 1,i,j,k,symi,symj,symk,rc1,                                      &
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,                     &
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,                     &
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
!
       if (typt3.eq.3) then
!4.4.2add part T2 . U
       call t3sgl(Work(iOff),wrksize,                                   &
     & mapdv,symijk,mapdfk4,mapifk4,mapdfk3,mapifk3,                    &
     & mapdt22,mapit22,mapdt21,mapit21,                                 &
     & 1,i,j,k,symi,symj,symk,rc1,                                      &
     & mapdm1,mapim1,possm10,mapdh1,mapih1,possh10,                     &
     & mapdm2,mapim2,possm20,mapdh2,mapih2,possh20,                     &
     & mapdm3,mapim3,possm30,mapdh3,mapih3,possh30)
       end if
!
!
!4.5  divide by denominators and calc energy contribution
!
!
!4.5.1divide by den.
       call t3div(Work(iOff),wrksize,                                   &
     & mapdw,mapdv,symijk,mapddp2,mapidp2,mapddp1,                      &
     & mapidp1,1,i,j,k,symi,symj,symk,ec,rc1)
!
!4.5.2add energy contribution
       ebbb=ebbb+ec
!
       end if
!
!
 1000   continue
 1100   continue
!
!par   Separate printing of partial energies e... are
!      useful only in serial run. For parallel run also
!      cycle over k is segmented (via paralelization),
!      so e... are not complete contributions (only sum
!      over all nodes have some sense). Thus, these values
!      in parallel run are too dangerous to use separately
!      so their printout is supressed.
!
       if (nProcs.eq.1) then
         call t3wresult (symi,symj,i,j,eaaa(1),eaab(1),eabb(1),ebbb(1))
         if (fullprint.gt.1) then
           write(6,*) ' Eaaa =',eaaa
           write(6,*) ' Eaab =',eaab
           write(6,*) ' Eabb =',eabb
           write(6,*) ' Ebbb =',ebbb
         end if
       end if
!endpar
!
 1201   continue
 1200   continue
 1300   continue
 1400   continue
!
!
!o    ***** final section *****
!
!o.*        allreduced energy components
        call gadgop (eaaa,1,'+')
        call gadgop (eaab,1,'+')
        call gadgop (eabb,1,'+')
        call gadgop (ebbb,1,'+')
!stare  call MPI_ALLREDUCE (ebbb,ec,1,
!    c  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)
!       ebbb=ec
!
!
!o.*  type results
!
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
! Export a method and energy to the MOLCAS runfile
      Call Put_cArray('Relax Method','CCSDT   ',8)
      Call Store_Energies(1,eccsd+eaaa+eaab+eabb+ebbb,1)
!     Releasing the memory
      Call GetMem('CCT3','Free','Real',iOff,wrksize)
!
!
      ireturn=0
      return
      end
