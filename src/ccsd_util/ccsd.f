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
       subroutine ccsd(ireturn,run_triples)
c
c     program for CCSD
c
       use Para_Info, only: MyRank, nProcs
       Logical run_triples
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "filemgr.fh"
#include "parallel.fh"

#include "SysDef.fh"
c
c     work file declaration
       integer wrksize
#include "WrkSpc.fh"
c
c     variables
c
       integer rc,posst,niter,idum(1)
       integer lunabij1,lunabij2,lunabij3
       integer lunt2o1,lunt2o2,lunt2o3
       integer lunw3aaaa,lunw3baab,lunw3bbaa,lunw3bbbb,lunw3abba,
     & lunw3aabb
       integer keyext,keyexc,lunrst
       real*8 scalar,energy,energyold,diff,dum(1)
       integer diispointt(1:4)
       integer diispointr(1:4)
       integer lenv,lenn,possabstack,nabstack,nfree,inv4,infree

       real*8 E1aa,E1bb,E2aaaa,E2bbbb,E2abab
       real*8 pz1aa,pz1bb,pz2abab
CLD    real*8 pz1aa,pz1bb,pz2aaaa,pz2bbbb,pz2abab
c
c      parameters for par
       real*8 timtotcpu,timtotcpun,timtotwc,timtotwcn,timdifwc
       real*8 timtotit
       integer i
c
       call CWTime (timtotcpu,timtotwc)
       timtotit=timtotwc
c
       energy = 9999999
       energyold=99999
       fullprint=0
       If (iPrintLevel(-1).LE.0) fullprint=-1

c
cI    **************  start   section **************
c
cI.1.1read input data form INPDAT (reorg) and input file
c
       call reainput
c
cI.1.2write head to output file
c
       if (fullprint.ge.0) then
         call wrhead
         if (noccsd.eq.1) then
           write(6,'(6X,A)') 'NOSD key is turned on'
           write(6,'(6X,A)') 'CCSD step is skipped'
           goto 999
         end if
       write (6,*) ' nProcs, myRank',nProcs,myRank
       end if
c
cI.1.3check if there is something to do in CCSD
       call ccsd_exc (keyexc)
       if (keyexc.eq.0) then
         write(6,*) ' No determinants in CCSD expansion '
         write(6,*) ' Zero correlation energy           '
         write(6,*) ' CCSD step excluded                '
         if (run_triples) then
           write(6,*) ' Triples step excluded'
           run_triples=.false.
         end if
         goto 999
       else if (keyexc.eq.1) then
         write(6,*) ' Only monoexcitations in CC expansion '
         write(6,*) ' CCSD step excluded                '
         if (run_triples) then
           write(6,*) ' Triples step excluded'
           run_triples=.false.
         end if
         goto 999
       end if
c
cI.1.4def parametrs for filemgr
c
       call mkfilemgrcom
c
cI.2.1 calc. work space requirements to fix and help mediates
c
       call initfiles (wrksize,lenv,lenn)
       if (fullprint.ge.0)
     &     write(6,*) ' Basic Work space requirements    :',wrksize
c
cI.2.3 allocate work space
c
cI.2.3.1  check free space
       Call GetMem('CCSD','Max','Real',maxspace,maxspace)
       maxspace=maxspace-8
       if (fullprint.ge.0)
     &     write (6,*) ' Max Size',maxspace
c
       if (maxspace.lt.wrksize) then
         write(6,*) ' Allocation of work space failed'
         write(6,*) ' Increase the value of the variable MOLCAS_MEM'
         Call Abend()
       end if
c
cI.2.3.2  fix, where AB stack will be located
       nfree=maxspace-wrksize
       inv4=int(lenv/lenn)
       infree=int(nfree/lenn)
       if (infree.gt.int(inv4/2)) then
c      there is enough free room for AB_stack
         possabstack=wrksize+1
         if (inv4.ge.infree) then
         nabstack=infree
         else
         nabstack=inv4
         end if
         wrksize=wrksize+nabstack*lenn
       else
c      there is not enough free room for AB_stack, will be build in
c      2.nd half of V4
         nabstack=int(inv4/2)
         possabstack=possv40+nabstack*lenn
       end if
c
       if (fullprint.ge.2) then
       write(6,*) ' Dimension of AB stack            :',nabstack
       end if
       if (fullprint.ge.0)
     &     write(6,*) ' Final Work space requirements    :',wrksize
c
cI.2.3.3  check free space
       Call GetMem('CCSD','Allo','Real',iOff,wrksize)
c
cI.2.4 set wrk = 0
       call mv0zero (wrksize,wrksize,Work(iOff))
       if (fullprint.ge.0)
     &     write(6,*) ' Allocation of work space   : Done'
c
cI.2.5read static integrals from INTSTA (reorg) file
c
       call reaintsta (Work(iOff),wrksize)
c
cI.2.6divide fok to faa,fai,fii and dp
c
       call divfok (Work(iOff),wrksize,
     & mapdn,mapin,possn0,mapdp,mapip,possp0,
     & mapdfk1,mapifk1,possfk10,mapdfk2,mapifk2,possfk20,
     & mapdfk3,mapifk3,possfk30,mapdfk4,mapifk4,possfk40,
     & mapdfk5,mapifk5,possfk50,mapdfk6,mapifk6,possfk60,
     & mapddp1,mapidp1,possdp10,mapddp2,mapidp2,possdp20,rc)
c
cI.2.7open 3 files for <ab|ij>aaaa,<ab|ij>bbbb,<ab|ij>abab
       call filemanager (1,lunabij1,rc)
       call filemanager (1,lunabij2,rc)
       call filemanager (1,lunabij3,rc)
c
cI.2.8open 3 files for T2oaaaa,T2obbbb,T2oabab
       call filemanager (1,lunt2o1,rc)
       call filemanager (1,lunt2o2,rc)
       call filemanager (1,lunt2o3,rc)
c
       if (keyrst.ne.0) then
cI.2.91        open fil for restart informations - filerst
       lunrst=16
       call filemanager (4,lunrst,rc)
       end if
c
       if (keyrst.eq.2) then
cI.2.92        read informations from restart file
c
c     get T1aa
       call getmediate (Work(iOff),wrksize,
     & lunrst,posst110,mapdt11,mapit11,rc)
c     get T1bb
       call getmediate (Work(iOff),wrksize,
     & lunrst,posst120,mapdt12,mapit12,rc)
c
c     get and write T2aaaa
       call getmediate (Work(iOff),wrksize,
     & lunrst,possv40,mapdv4,mapiv4,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunt2o1,mapdv4,mapiv4,rc)
c     get and write T2bbbb
       call getmediate (Work(iOff),wrksize,
     & lunrst,possv40,mapdv4,mapiv4,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunt2o2,mapdv4,mapiv4,rc)
c     get and write T2abab
       call getmediate (Work(iOff),wrksize,
     & lunrst,possv40,mapdv4,mapiv4,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunt2o3,mapdv4,mapiv4,rc)
c
       if (iokey.eq.1) then
c      Fortran IO
       read (lunrst,end=97) energyold,niter
c
       else
c      MOLCAS IO
       call ddafile (lunrst,2,dum,1,daddr(lunrst))
       energyold=dum(1)
       call idafile (lunrst,2,idum,1,daddr(lunrst))
       niter=idum(1)
       end if
       goto 98
c
c     case when energy and niter is not in save file
CFUE 97     energyold=0.0d0
 97     energyold=Escf
       niter=0
       write(6,*) ' ENERGY AND NIT WAS NOT IN SAVE FILE, CHANGED TO 0'
c
 98     niter=niter+1
c
       else
cI.2.93 T2naaaa,T2nbbbb and T2nabab are 0, we can write them as T2old
       call wrtmediate (Work(iOff),wrksize,
     & lunt2o1,mapdt21,mapit21,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunt2o2,mapdt22,mapit22,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunt2o3,mapdt23,mapit23,rc)
       end if
c
c
cI.2.10        write V1 (<ab||ij>aaaa) to lunabij1
c     V2 (<ab||ij>bbbb) to lunabij2
c     V3 (<ab||ij>abab) to lunabij3
c
       call wrtmediate (Work(iOff),wrksize,
     & lunabij1,mapdv1,mapiv1,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunabij2,mapdv2,mapiv2,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunabij3,mapdv3,mapiv3,rc)
c
cI.2.11        open lunext files if extrapolation is used
c
       if (yesext.eq.1) then
c     keyext=firstext
c
c*    open files for T and R
       call diisof (diispointt,cycext)
       call diisof (diispointr,cycext)
c
c*    open file for E (extrapolate amplitudes)
       call filemanager (1,lune,rc)
c
c*    write T21,T22,T23,T11,T12 to E (they are zero or restart ones)
c     T2aaaa
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt21,mapit21,rc)
c     T2bbbb
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt22,mapit22,rc)
c     T2abab
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt23,mapit23,rc)
c     T1aa
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt11,mapit11,rc)
c     T1bb
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt12,mapit12,rc)
c
c*
       keyext=0
       end if
c
c*     write head for short + medium printing
       if (fullprint.ge.0 .and. fullprint.lt.2) then
         write(6,*)
         write(6,'(6X,A)')
     &   'Iteration       Total enegy      Corr. energy      Difference'
       end if
c
c
cII   **************  Iteration cycle   **************
c
       if (keyrst.ne.2) then
         niter=1
       end if
c
       call distnodes
c
c
 1     if ((fullprint.ge.1).and.(myRank.eq.0).and.(nProcs.gt.1)) then
         write (6,*) 'Npocab:',nprocab
         do i=1,nprocab
           write (6,*) 'Proc  :',i,idab(i),ideffab(i)
         end do
       write (6,*) 'IDaaaa:',idaaaa
       write (6,*) 'IDaabb:',idaabb
       write (6,*) 'IDabba:',idabba
       write (6,*) 'IDbbbb:',idbbbb
       write (6,*) 'IDbbaa:',idbbaa
       write (6,*) 'IDbaab:',idbaab
       write (6,*) 'IDfin :',idfin
       end if
cpar
c    set ididle,idtmab =0
       do i=1,nProcs
        ididle(i)=0.0d0
        idtmab(i)=0.0d0
       end do
cend par
c
       call init (Work(iOff),wrksize,
     & lunabij1,lunabij2,lunabij3)
       call CWTime (timtotcpun,timtotwcn)
ctmp   timdifcpu=timtotcpun-timtotcpu
       timdifwc=timtotwcn-timtotwc
       timtotcpu=timtotcpun
       timtotwc=timtotwcn
       if (fullprint.ge.1) then
         write(6,*) ' Initial section done',myRank,timdifwc,
     &                                     timtotwcn-timtotit
       timtotit=timtotwcn
       end if
c
       call sumoverab (Work(iOff),wrksize,
     & lunt2o1,lunt2o2,lunt2o3,nabstack,possabstack,niter)
       call CWTime (timtotcpun,timtotwcn)
ctmp   timdifcpu=timtotcpun-timtotcpu
       timdifwc=timtotwcn-timtotwc
       timtotcpu=timtotcpun
       timtotwc=timtotwcn
       if (fullprint.ge.1) then
         write(6,*) ' Sumation over ab done',myRank,timdifwc
       end if
c
cpar   store time spend in sumoverab process for this node
       idtmab(myRank+1)=timdifwc
c
       call sumovera (Work(iOff),wrksize,
     & lunt2o1,lunt2o2,lunt2o3,
     & lunw3aaaa,lunw3baab,lunw3bbaa,lunw3bbbb,lunw3abba,lunw3aabb)
       call CWTime (timtotcpun,timtotwcn)
ctmp   timdifcpu=timtotcpun-timtotcpu
       timdifwc=timtotwcn-timtotwc
       timtotcpu=timtotcpun
       timtotwc=timtotwcn
       if (fullprint.ge.1) then
         write(6,*) ' Sumation over a done',myRank,timdifwc
       end if
c
       call intermezzo (Work(iOff),wrksize,
     & lunw3aaaa,lunw3bbbb,lunw3abba,
     & lunw3baab,lunw3aabb,lunw3bbaa,lunt2o1,lunt2o2,lunt2o3,
     & lunabij1,lunabij2,lunabij3)
       call CWTime (timtotcpun,timtotwcn)
ctmp   timdifcpu=timtotcpun-timtotcpu
       timdifwc=timtotwcn-timtotwc
       timtotcpu=timtotcpun
       timtotwc=timtotwcn
       if (fullprint.ge.1) then
         write(6,*) ' Internal section done',myRank,timdifwc
       end if
c
       call finale (Work(iOff),wrksize,
     & lunabij1,lunabij2,lunabij3,
     & lunt2o1,lunt2o2,lunt2o3)
       call CWTime (timtotcpun,timtotwcn)
ctmp   timdifcpu=timtotcpun-timtotcpu
       timdifwc=timtotwcn-timtotwc
       timtotcpu=timtotcpun
       timtotwc=timtotwcn
       if (fullprint.ge.1) then
         write(6,*) ' Final section done',myRank,timdifwc
       end if
c
#ifdef _MOLCAS_MPP_
cpar
       call joinamplitudes (Work(iOff),wrksize)
       call CWTime (timtotcpun,timtotwcn)
ctmp   timdifcpu=timtotcpun-timtotcpu
       timdifwc=timtotwcn-timtotwc
       timtotcpu=timtotcpun
       timtotwc=timtotwcn
       if (fullprint.ge.1) then
         write(6,*) ' Amplitudes allreduced',myRank,timdifwc
       end if
c
c      store differntial wc time at this point to calculate idle time
c      and redefine ideffab for next step
       ididle(myRank+1)=timdifwc
       if (niter.gt.0) then
         call redef
         if ((myRank.eq.0).and.(fullprint.ge.1)) then
         write (6,997) ' ABTim',(idtmab(i),i=1,nProcs)
           write (6,997) ' IDtim',(ididle(i),i=1,nProcs)
           write (6,997) ' EFF  ',(ideffab(i),i=1,nprocab)
997        format (a6,2x,16(f7.2,1x))
       end if
       end if
cendpar
#endif
c
cIII  *********** division, SA, extrapolation, energy, tests *************
c
c
c3.1  division by denominators (divided amplitudes are stored in the same files
c     like T1n and T2n
c
c3.1.1div. T1n
       call divt (Work(iOff),wrksize,
     & 2,mapdt13,mapit13,mapddp1,mapidp1,mapddp2,mapidp2,rc)
       call divt (Work(iOff),wrksize,
     & 2,mapdt14,mapit14,mapddp1,mapidp1,mapddp2,mapidp2,rc)
c
c3.1.2div. T2n
       call divt (Work(iOff),wrksize,
     & 4,mapdt21,mapit21,mapddp1,mapidp1,mapddp2,mapidp2,rc)
       call divt (Work(iOff),wrksize,
     & 4,mapdt22,mapit22,mapddp1,mapidp1,mapddp2,mapidp2,rc)
       call divt (Work(iOff),wrksize,
     & 4,mapdt23,mapit23,mapddp1,mapidp1,mapddp2,mapidp2,rc)
c
c
c3.2  Spin adaptation
c
       if (keysa.gt.0) then
       call saamp (Work(iOff),wrksize,
     & keysa)
       end if
c
c
c3.3  extrapolation
       if (yesext.eq.1) then
c
c*    write new Tn (T21,T22,T23,T13,T14) to T stack
       call diiswa1 (Work(iOff),wrksize,
     & diispointt)
c
c*    calc Tn = Tn - E
       call calcr (Work(iOff),wrksize,
     & lune)
c
c*    write Tn=Tn-E to R stack
       call diiswa1 (Work(iOff),wrksize,
     & diispointr)
c
c*    do diis Tn = DIIS (Tprev) if necc.
       call diis (Work(iOff),wrksize,
     & diispointt,diispointr,keyext)
c
c*    write Tn = T(DIIS) to E
c     rewind lune
       call filemanager (2,lune,rc)
c     T2aaaa
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt21,mapit21,rc)
c     T2bbbb
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt22,mapit22,rc)
c     T2abab
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt23,mapit23,rc)
c     T1aa
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt13,mapit13,rc)
c     T1bb
       call wrtmediate (Work(iOff),wrksize,
     & lune,mapdt14,mapit14,rc)
c
c
       end if
c
c
c3.4  save restart informations - amplitudes
       if (keyrst.ne.0) then
       call saverest1 (Work(iOff),wrksize,
     & lunrst)
       end if
c
c
c3.5  write Tn into place of to
c
c3.5.1put t1naa -> t1oaa
       call map (Work(iOff),wrksize,
     & 2,1,2,0,0,mapdt13,mapit13,1,mapdt11,mapit11,posst110,
     &           posst,rc)
c
c3.5.2put t1nbb -> t1obb
       call map (Work(iOff),wrksize,
     & 2,1,2,0,0,mapdt14,mapit14,1,mapdt12,mapit12,posst120,
     &           posst,rc)
c
c3.5.3rewind lunt2o1 and write t2naaaa there
       call filemanager (2,lunt2o1,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunt2o1,mapdt21,mapit21,rc)
c
c3.5.4rewind lunt2o2 and write t2nbbbb there
       call filemanager (2,lunt2o2,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunt2o2,mapdt22,mapit22,rc)
c
c3.5.5rewind lunt2o3 and write t2nabab there
       call filemanager (2,lunt2o3,rc)
       call wrtmediate (Work(iOff),wrksize,
     & lunt2o3,mapdt23,mapit23,rc)
c
c     test zero
       call percentzero (Work(iOff),wrksize,
     & mapdt13,pz1aa)
       call percentzero (Work(iOff),wrksize,
     & mapdt14,pz1bb)
c      call percentzero (Work(iOff),wrksize,
c    & mapdt21,pz2aaaa)
c      call percentzero (Work(iOff),wrksize,
c    & mapdt22,pz2bbbb)
       call percentzero (Work(iOff),wrksize,
     & mapdt23,pz2abab)
c
c
c3.6  calc energy
CFUE   energy=0.0d0
       energy=Escf
       scalar=0.0d0
c
c3.6.1calc fai(a,i)aa . t1(a,i)aa
       call multdot (Work(iOff),wrksize,
     & 2,mapdfk3,mapifk3,1,mapdt13,mapit13,1,scalar,rc)
       energy=energy+scalar
       E1aa = scalar
CFUE   write(6,22) scalar
c
c3.6.2calc fai(a,i)bb . t1(a,i)bb
       call multdot (Work(iOff),wrksize,
     & 2,mapdfk4,mapifk4,1,mapdt14,mapit14,1,scalar,rc)
       energy=energy+scalar
       E1bb = scalar
CFUE   write(6,22) scalar
c
c3.6.3get <ab||ij>aaaa into V1 and calc V1(abij) . Tau(abij)aaaa
       call filemanager (2,lunabij1,rc)
       call getmediate (Work(iOff),wrksize,
     & lunabij1,possv10,mapdv1,mapiv1,rc)
       call mktau (Work(iOff),wrksize,
     & mapdt21,mapit21,mapdt13,mapit13,mapdt13,mapit13,
     &             1.0d0,rc)
       call multdot (Work(iOff),wrksize,
     & 4,mapdv1,mapiv1,1,mapdt21,mapit21,1,scalar,rc)
       energy=energy+scalar
       E2aaaa = scalar
CFUE   write(6,22) scalar
c
c3.6.4get <ab||ij>bbbb into V1 and calc V1(abij) . Tau(abij)bbbb
       call filemanager (2,lunabij2,rc)
       call getmediate (Work(iOff),wrksize,
     & lunabij2,possv10,mapdv1,mapiv1,rc)
       call mktau (Work(iOff),wrksize,
     & mapdt22,mapit22,mapdt14,mapit14,mapdt14,mapit14,
     &             1.0d0,rc)
       call multdot (Work(iOff),wrksize,
     & 4,mapdv1,mapiv1,1,mapdt22,mapit22,1,scalar,rc)
       energy=energy+scalar
       E2bbbb = scalar
CFUE   write(6,22) scalar
c
c3.6.5get <ab||ij>abab into V1 and calc V1(abij) . Tau(abij)abab
       call filemanager (2,lunabij3,rc)
       call getmediate (Work(iOff),wrksize,
     & lunabij3,possv10,mapdv1,mapiv1,rc)
       call mktau (Work(iOff),wrksize,
     & mapdt23,mapit23,mapdt13,mapit13,mapdt14,mapit14,
     &             1.0d0,rc)
       call multdot (Work(iOff),wrksize,
     & 4,mapdv1,mapiv1,1,mapdt23,mapit23,1,scalar,rc)
       energy=energy+scalar
       E2abab = scalar
CFUE   write(6,22) scalar
CFUE   write(6,22) energy
C22    format (2x,f20.15)
c
c
c3.7  save restart informations - energy, iteration cycle
       if (keyrst.ne.0) then
       call saverest2 (lunrst,energy,niter,iokey,daddr(lunrst))
       end if
c
ctmp   call timing (timtotcpu,timdifcpu,timtotwc,timdifwc)
       call CWTime (timtotcpun,timtotwcn)
ctmp   timdifcpu=timtotcpun-timtotcpu
       timdifwc=timtotwcn-timtotwc
       timtotcpu=timtotcpun
       timtotwc=timtotwcn
       if (fullprint.ge.1) then
         write(6,*) ' Adapt., Extrap. and Energy evaluation done ',
     c               myRank,timdifwc,timtotwc
       end if
c
c
c3.8  test of convergence or termination
c
       diff=energy-energyold
       If (niter.eq.1) diff=E1aa+E1bb+E2aaaa+E2bbbb+E2abab
       energyold=energy
c
       if (fullprint.ge.0 .and. fullprint.le.1) then
c      reduced printing
         write(6,'(6X,I5,5X,F17.8,2X,2(F15.8,2X))')
     &   niter,energy,E1aa+E1bb+E2aaaa+E2bbbb+E2abab,diff
c
       else if (fullprint.ge.2) then
c      full printing
         write(6,*)
         write(6,'(6X,A,I4)') 'Iteration No        :',niter
         write(6,'(6X,A,2F24.13)') 'Total energy (diff) :',energy,diff
         write(6,'(6X,A,F24.13)') 'Correlation energy  :',
     &                         E1aa+E1bb+E2aaaa+E2bbbb+E2abab
         write(6,'(6X,A,F24.13)') 'Reference energy    :',Escf
         write(6,'(6X,A,F17.8)')  'E1aa   contribution :',E1aa
         write(6,'(6X,A,F17.8)')  'E1bb   contribution :',E1bb
         write(6,'(6X,A,F17.8)')  'E2aaaa contribution :',E2aaaa
         write(6,'(6X,A,F17.8)')  'E2bbbb contribution :',E2bbbb
         write(6,'(6X,A,F17.8)')  'E2abab contribution :',E2abab
         write(6,'(6X,A)')
     &   '% of small amplitudes in  T1aa   T1bb   T2aaaa T2bbbb T2abab'
         write(6,'(31X,5(F5.1,2X))')
     &        pz1aa,pz1bb,pz2abab
c
       end if
       Call Add_Info('E_CCSD',[Energy],1,8)
c
       niter=niter+1
c
*--- the original code depends on print level!!!!!!!
*      if (abs(diff).le.ccconv .and. fullprint.ge.0) then
*        write(6,*) '     Convergence after ',niter,
*    &   ' Iterations'
*      else if (niter.gt.maxiter .and. fullprint.ge.0) then
*        write(6,*) '     Convergence not reached'
*      else if (niter.le.maxiter)then
*        goto 1
*      end if
       if (abs(diff).le.ccconv) then
         write(6,*) '     Convergence after ',niter,
     &   ' Iterations'
       else if (niter.gt.maxiter) then
         write(6,*) '     Convergence not reached'
       else if (niter.le.maxiter)then
         goto 1
       end if
       IF (fullprint.ge.0) THEN
         write(6,*)
         write(6,*)
c
c
c**   **************     Finitto      **************
c
         write(6,'(6X,A,2F17.8)') 'Total energy (diff) :',energy,diff
         write(6,'(6X,A,F24.13)') 'Correlation energy  :',
     &                         E1aa+E1bb+E2aaaa+E2bbbb+E2abab
         write(6,'(6X,A,F24.13)') 'Reference energy    :',Escf
         write(6,'(6X,A,F17.8)')  'E1aa   contribution :',E1aa
         write(6,'(6X,A,F17.8)')  'E1bb   contribution :',E1bb
         write(6,'(6X,A,F17.8)')  'E2aaaa contribution :',E2aaaa
         write(6,'(6X,A,F17.8)')  'E2bbbb contribution :',E2bbbb
         write(6,'(6X,A,F17.8)')  'E2abab contribution :',E2abab
         write(6,*)
         write(6,*)
       ENDIF
c Export a method and energy to the MOLCAS runfile
       Call Put_cArray('Relax Method','CCSDT   ',8)
       Call Store_Energies(1,[Energy],1)
c
c4.0  type 5 maximal elements + euclidian norms in each type of amplitudes
c
c4.0.1T1aa
       call max5 (Work(iOff),wrksize,
     & 2,mapdt13,mapit13,'T1aa    ')
c
c4.0.2T1bb
       call max5 (Work(iOff),wrksize,
     & 2,mapdt14,mapit14,'T1bb    ')
c
c4.0.3T2aaaa
       call max5 (Work(iOff),wrksize,
     & 4,mapdt21,mapit21,'T2aaaa  ')
c
c4.0.4T2bbbb
       call max5 (Work(iOff),wrksize,
     & 4,mapdt22,mapit22,'T2bbbb  ')
c
c4.0.5T2abab
       call max5 (Work(iOff),wrksize,
     & 4,mapdt23,mapit23,'T2abab  ')
c
c4.1  close lunabij1,2,3
       call filemanager (3,lunabij1,rc)
       call filemanager (3,lunabij2,rc)
       call filemanager (3,lunabij3,rc)
c
c4.2  close lunt2o1,2,3
       call filemanager (3,lunt2o1,rc)
       call filemanager (3,lunt2o2,rc)
       call filemanager (3,lunt2o3,rc)
c
c4.3  close lunext files if extrapolation is used
       if (yesext.eq.1) then
       call diiscf (diispointt,cycext)
       call diiscf (diispointr,cycext)
       call filemanager (3,lune,rc)
       end if
c
c4.4  close lunrst file if restart informations was stored
       if (keyrst.ne.0) then
       call filemanager (3,lunrst,rc)
       end if
c global ccsd synchronization point for parallel runs
       Call GASync()

c      Releasing the memory
       Call GetMem('CCSD','Free','Real',iOff,wrksize)
c
999    If (fullprint.ge.0) then
         write(6,*)
         write(6,'(6X,A)') 'Happy Landing!'
         write(6,*)
       EndIf
c
       ireturn=0
       return
       end
c
c     ---------------
c
       subroutine percentzero (wrk,wrksize,
     & mapd,pz)
c
c     this routine test % of small elements in meditate, decribed by mpd
c
c     mapd - direct map of required mediate (I)
c
#include "wrk.fh"
       integer mapd(0:512,1:6)
       real*8  pz
c
c     help variables
c
       integer poss,length
       integer nhelp,nzero
       real*8 zerolim
c
c     def length, poss, zerolim
c
       poss=mapd(1,1)
       nhelp=mapd(0,5)
       length=mapd(nhelp,1)+mapd(nhelp,2)-mapd(1,1)
       zerolim=1.0d-6
c
       if (length.gt.0) then
       nzero=0
       do 100 nhelp=poss,poss+length-1
       if (abs(wrk(nhelp)).lt.zerolim) then
       nzero=nzero+1
       end if
 100    continue
       pz = dble(100*nzero)/dble(length)
       else
       pz=1.0d0
       end if
c
       return
       end
c
c     ---------------
c
        subroutine ccsd_exc (key)
c
c       check, if there is atleast one determinant in CCSD expansion
c       key=0 - no determinant in expansion
c           1 - only monoexcitations in expansion
c           2 - both mono and biexcitations in expansion
c
        implicit none
#include "ccsd1.fh"
        integer key
c
c       help variables
        integer isym,jsym,ijsym,asym,bsym,nij,nab
        integer naa,nbb,naaaa,nbbbb,nabab
c
c
c1.1    calc # of monoexcitations
c       taking into account also symmetry

        naa=0
        nbb=0
        do isym=1,nsym
          asym=isym
          naa=naa+noa(isym)*nva(asym)
          nbb=nbb+nob(isym)*nvb(asym)
        end do
c
c1.2    calc # of biexcitation
c       taking into account also symmetry
c
        naaaa=0
        do isym=1,nsym
        do jsym=1,isym
        ijsym=mmul(isym,jsym)
        if (isym.eq.jsym) then
          nij=noa(isym)*(noa(isym)-1)/2
        else
          nij=noa(isym)*noa(jsym)
        end if
          do asym=1,nsym
          bsym=mmul(ijsym,asym)
          if (bsym.lt.asym) then
            nab=nva(asym)*nva(bsym)
          else if (bsym.eq.asym) then
            nab=nva(asym)*(nva(asym)-1)/2
          else
            nab=0
          end if
          naaaa=naaaa+nij*nab
          end do
        end do
        end do
c
        nbbbb=0
        do isym=1,nsym
        do jsym=1,isym
        ijsym=mmul(isym,jsym)
        if (isym.eq.jsym) then
          nij=nob(isym)*(nob(isym)-1)/2
        else
          nij=nob(isym)*nob(jsym)
        end if
          do asym=1,nsym
          bsym=mmul(ijsym,asym)
          if (bsym.lt.asym) then
            nab=nvb(asym)*nvb(bsym)
          else if (bsym.eq.asym) then
            nab=nvb(asym)*(nvb(asym)-1)/2
          else
            nab=0
          end if
          nbbbb=nbbbb+nij*nab
          end do
        end do
        end do
c
        nabab=0
        do isym=1,nsym
        do jsym=1,isym
        ijsym=mmul(isym,jsym)
        nij=noa(isym)*nob(jsym)
          do asym=1,nsym
          bsym=mmul(ijsym,asym)
          nab=nva(asym)*nvb(bsym)
          nabab=nabab+nij*nab
          end do
        end do
        end do
c
c
c2      set key
c
        if ((naaaa+nbbbb+nabab).eq.0) then
          if ((naa+nbb).eq.0) then
            key=0
          else
            key=1
          end if
        else
          key=2
        end if
c
c
        return
        end
