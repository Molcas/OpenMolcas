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

subroutine ccsd(ireturn,run_triples)
! program for CCSD

use Para_Info, only: MyRank, nProcs
logical run_triples
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "filemgr.fh"
#include "parallel.fh"
#include "SysDef.fh"
! work file declaration
integer wrksize
#include "WrkSpc.fh"
! variables
integer rc, posst, niter, idum(1)
integer lunabij1, lunabij2, lunabij3
integer lunt2o1, lunt2o2, lunt2o3
integer lunw3aaaa, lunw3baab, lunw3bbaa, lunw3bbbb, lunw3abba, lunw3aabb
integer keyext, keyexc, lunrst
real*8 scalar, energy, energyold, diff, dum(1)
integer diispointt(1:4)
integer diispointr(1:4)
integer lenv, lenn, possabstack, nabstack, nfree, inv4, infree
real*8 E1aa, E1bb, E2aaaa, E2bbbb, E2abab
real*8 pz1aa, pz1bb, pz2abab
!LD real*8 pz1aa,pz1bb,pz2aaaa,pz2bbbb,pz2abab
! parameters for par
real*8 timtotcpu, timtotcpun, timtotwc, timtotwcn, timdifwc
real*8 timtotit
integer i

call CWTime(timtotcpu,timtotwc)
timtotit = timtotwc

energy = 9999999
energyold = 99999
fullprint = 0
if (iPrintLevel(-1) <= 0) fullprint = -1

!I **************  start   section **************

!I.1.1 read input data form INPDAT (reorg) and input file

call reainput()

!I.1.2 write head to output file

if (fullprint >= 0) then
  call wrhead()
  if (noccsd == 1) then
    write(6,'(6X,A)') 'NOSD key is turned on'
    write(6,'(6X,A)') 'CCSD step is skipped'
    goto 999
  end if
  write(6,*) ' nProcs, myRank',nProcs,myRank
end if

!I.1.3 check if there is something to do in CCSD
call ccsd_exc(keyexc)
if (keyexc == 0) then
  write(6,*) ' No determinants in CCSD expansion '
  write(6,*) ' Zero correlation energy           '
  write(6,*) ' CCSD step excluded                '
  if (run_triples) then
    write(6,*) ' Triples step excluded'
    run_triples = .false.
  end if
  goto 999
else if (keyexc == 1) then
  write(6,*) ' Only monoexcitations in CC expansion '
  write(6,*) ' CCSD step excluded                '
  if (run_triples) then
    write(6,*) ' Triples step excluded'
    run_triples = .false.
  end if
  goto 999
end if

!I.1.4 def parametrs for filemgr

call mkfilemgrcom()

!I.2.1 calc. work space requirements to fix and help mediates

call initfiles(wrksize,lenv,lenn)
if (fullprint >= 0) write(6,*) ' Basic Work space requirements    :',wrksize

!I.2.3 allocate work space

!I.2.3.1 check free space
call GetMem('CCSD','Max','Real',idum(1),maxspace)
maxspace = maxspace-8
if (fullprint >= 0) write(6,*) ' Max Size',maxspace

if (maxspace < wrksize) then
  write(6,*) ' Allocation of work space failed'
  write(6,*) ' Increase the value of the variable MOLCAS_MEM'
  call Abend()
end if

!I.2.3.2  fix, where AB stack will be located
nfree = maxspace-wrksize
inv4 = int(lenv/lenn)
infree = int(nfree/lenn)
if (infree > int(inv4/2)) then
  ! there is enough free room for AB_stack
  possabstack = wrksize+1
  if (inv4 >= infree) then
    nabstack = infree
  else
    nabstack = inv4
  end if
  wrksize = wrksize+nabstack*lenn
else
  ! there is not enough free room for AB_stack, will be build in
  ! 2nd half of V4
  nabstack = int(inv4/2)
  possabstack = possv40+nabstack*lenn
end if

if (fullprint >= 2) then
  write(6,*) ' Dimension of AB stack            :',nabstack
end if
if (fullprint >= 0) write(6,*) ' Final Work space requirements    :',wrksize

!I.2.3.3 check free space
call GetMem('CCSD','Allo','Real',iOff,wrksize)

!I.2.4 set wrk = 0
call mv0zero(wrksize,wrksize,Work(iOff))
if (fullprint >= 0) write(6,*) ' Allocation of work space   : Done'

!I.2.5 read static integrals from INTSTA (reorg) file

call reaintsta(Work(iOff),wrksize)

!I.2.6 divide fok to faa,fai,fii and dp

call divfok(Work(iOff),wrksize,mapdn,mapin,possn0,mapdp,mapip,possp0,mapdfk1,mapifk1,possfk10,mapdfk2,mapifk2,possfk20,mapdfk3, &
            mapifk3,possfk30,mapdfk4,mapifk4,possfk40,mapdfk5,mapifk5,possfk50,mapdfk6,mapifk6,possfk60,mapddp1,mapidp1,possdp10, &
            mapddp2,mapidp2,possdp20,rc)

!I.2.7 open 3 files for <ab|ij>aaaa,<ab|ij>bbbb,<ab|ij>abab
call filemanager(1,lunabij1,rc)
call filemanager(1,lunabij2,rc)
call filemanager(1,lunabij3,rc)

!I.2.8 open 3 files for T2oaaaa,T2obbbb,T2oabab
call filemanager(1,lunt2o1,rc)
call filemanager(1,lunt2o2,rc)
call filemanager(1,lunt2o3,rc)

if (keyrst /= 0) then
  !I.2.9.1 open fil for restart informations - filerst
  lunrst = 16
  call filemanager(4,lunrst,rc)
end if

if (keyrst == 2) then
  !I.2.9.2 read informations from restart file

  ! get T1aa
  call getmediate(Work(iOff),wrksize,lunrst,posst110,mapdt11,mapit11,rc)
  ! get T1bb
  call getmediate(Work(iOff),wrksize,lunrst,posst120,mapdt12,mapit12,rc)

  ! get and write T2aaaa
  call getmediate(Work(iOff),wrksize,lunrst,possv40,mapdv4,mapiv4,rc)
  call wrtmediate(Work(iOff),wrksize,lunt2o1,mapdv4,mapiv4,rc)
  ! get and write T2bbbb
  call getmediate(Work(iOff),wrksize,lunrst,possv40,mapdv4,mapiv4,rc)
  call wrtmediate(Work(iOff),wrksize,lunt2o2,mapdv4,mapiv4,rc)
  ! get and write T2abab
  call getmediate(Work(iOff),wrksize,lunrst,possv40,mapdv4,mapiv4,rc)
  call wrtmediate(Work(iOff),wrksize,lunt2o3,mapdv4,mapiv4,rc)

  if (iokey == 1) then
    ! Fortran IO
    read(lunrst,end=97) energyold,niter

  else
    ! MOLCAS IO
    call ddafile(lunrst,2,dum,1,daddr(lunrst))
    energyold = dum(1)
    call idafile(lunrst,2,idum,1,daddr(lunrst))
    niter = idum(1)
  end if
  goto 98

  ! case when energy and niter is not in save file
  !FUE 97 energyold = 0.0d0
  97 continue
  energyold = Escf
  niter = 0
  write(6,*) ' ENERGY AND NIT WAS NOT IN SAVE FILE, CHANGED TO 0'

  98 continue
  niter = niter+1

else
  !I.2.9.3 T2naaaa,T2nbbbb and T2nabab are 0, we can write them as T2old
  call wrtmediate(Work(iOff),wrksize,lunt2o1,mapdt21,mapit21,rc)
  call wrtmediate(Work(iOff),wrksize,lunt2o2,mapdt22,mapit22,rc)
  call wrtmediate(Work(iOff),wrksize,lunt2o3,mapdt23,mapit23,rc)
end if

!I.2.10 write V1 (<ab||ij>aaaa) to lunabij1
!     V2 (<ab||ij>bbbb) to lunabij2
!     V3 (<ab||ij>abab) to lunabij3

call wrtmediate(Work(iOff),wrksize,lunabij1,mapdv1,mapiv1,rc)
call wrtmediate(Work(iOff),wrksize,lunabij2,mapdv2,mapiv2,rc)
call wrtmediate(Work(iOff),wrksize,lunabij3,mapdv3,mapiv3,rc)

!I.2.11 open lunext files if extrapolation is used

if (yesext == 1) then
  !keyext = firstext

  !* open files for T and R
  call diisof(diispointt,cycext)
  call diisof(diispointr,cycext)

  !* open file for E (extrapolate amplitudes)
  call filemanager(1,lune,rc)

  !* write T21,T22,T23,T11,T12 to E (they are zero or restart ones)
  ! T2aaaa
  call wrtmediate(Work(iOff),wrksize,lune,mapdt21,mapit21,rc)
  ! T2bbbb
  call wrtmediate(Work(iOff),wrksize,lune,mapdt22,mapit22,rc)
  ! T2abab
  call wrtmediate(Work(iOff),wrksize,lune,mapdt23,mapit23,rc)
  ! T1aa
  call wrtmediate(Work(iOff),wrksize,lune,mapdt11,mapit11,rc)
  ! T1bb
  call wrtmediate(Work(iOff),wrksize,lune,mapdt12,mapit12,rc)

  keyext = 0
end if

!* write head for short + medium printing
if ((fullprint >= 0) .and. (fullprint < 2)) then
  write(6,*)
  write(6,'(6X,A)') 'Iteration       Total enegy      Corr. energy      Difference'
end if

!II **************  Iteration cycle   **************

if (keyrst /= 2) then
  niter = 1
end if

call distnodes()

1 continue
if ((fullprint >= 1) .and. (myRank == 0) .and. (nProcs > 1)) then
  write(6,*) 'Npocab:',nprocab
  do i=1,nprocab
    write(6,*) 'Proc  :',i,idab(i),ideffab(i)
  end do
  write(6,*) 'IDaaaa:',idaaaa
  write(6,*) 'IDaabb:',idaabb
  write(6,*) 'IDabba:',idabba
  write(6,*) 'IDbbbb:',idbbbb
  write(6,*) 'IDbbaa:',idbbaa
  write(6,*) 'IDbaab:',idbaab
  write(6,*) 'IDfin :',idfin
end if
!par
! set ididle,idtmab =0
do i=1,nProcs
  ididle(i) = 0.0d0
  idtmab(i) = 0.0d0
end do
!end par

call init(Work(iOff),wrksize,lunabij1,lunabij2,lunabij3)
call CWTime(timtotcpun,timtotwcn)
!tmp timdifcpu = timtotcpun-timtotcpu
timdifwc = timtotwcn-timtotwc
timtotcpu = timtotcpun
timtotwc = timtotwcn
if (fullprint >= 1) then
  write(6,*) ' Initial section done',myRank,timdifwc,timtotwcn-timtotit
  timtotit = timtotwcn
end if

call sumoverab(Work(iOff),wrksize,lunt2o1,lunt2o2,lunt2o3,nabstack,possabstack,niter)
call CWTime(timtotcpun,timtotwcn)
!tmp timdifcpu = timtotcpun-timtotcpu
timdifwc = timtotwcn-timtotwc
timtotcpu = timtotcpun
timtotwc = timtotwcn
if (fullprint >= 1) then
  write(6,*) ' Summation over ab done',myRank,timdifwc
end if

!par store time spend in sumoverab process for this node
idtmab(myRank+1) = timdifwc

call sumovera(Work(iOff),wrksize,lunt2o1,lunt2o2,lunt2o3,lunw3aaaa,lunw3baab,lunw3bbaa,lunw3bbbb,lunw3abba,lunw3aabb)
call CWTime(timtotcpun,timtotwcn)
!tmp timdifcpu = timtotcpun-timtotcpu
timdifwc = timtotwcn-timtotwc
timtotcpu = timtotcpun
timtotwc = timtotwcn
if (fullprint >= 1) then
  write(6,*) ' Summation over a done',myRank,timdifwc
end if

call intermezzo(Work(iOff),wrksize,lunw3aaaa,lunw3bbbb,lunw3abba,lunw3baab,lunw3aabb,lunw3bbaa,lunt2o1,lunt2o2,lunt2o3,lunabij1, &
                lunabij2,lunabij3)
call CWTime(timtotcpun,timtotwcn)
!tmp timdifcpu = timtotcpun-timtotcpu
timdifwc = timtotwcn-timtotwc
timtotcpu = timtotcpun
timtotwc = timtotwcn
if (fullprint >= 1) then
  write(6,*) ' Internal section done',myRank,timdifwc
end if

call finale(Work(iOff),wrksize,lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
call CWTime(timtotcpun,timtotwcn)
!tmp timdifcpu = timtotcpun-timtotcpu
timdifwc = timtotwcn-timtotwc
timtotcpu = timtotcpun
timtotwc = timtotwcn
if (fullprint >= 1) then
  write(6,*) ' Final section done',myRank,timdifwc
end if

#ifdef _MOLCAS_MPP_
!par
call joinamplitudes(Work(iOff),wrksize)
call CWTime(timtotcpun,timtotwcn)
!tmp timdifcpu = timtotcpun-timtotcpu
timdifwc = timtotwcn-timtotwc
timtotcpu = timtotcpun
timtotwc = timtotwcn
if (fullprint >= 1) then
  write(6,*) ' Amplitudes allreduced',myRank,timdifwc
end if

! store differential wc time at this point to calculate idle time
! and redefine ideffab for next step
ididle(myRank+1) = timdifwc
if (niter > 0) then
  call redef()
  if ((myRank == 0) .and. (fullprint >= 1)) then
    write(6,997) ' ABTim',(idtmab(i),i=1,nProcs)
    write(6,997) ' IDtim',(ididle(i),i=1,nProcs)
    write(6,997) ' EFF  ',(ideffab(i),i=1,nprocab)
997 format(a6,2x,16(f7.2,1x))
  end if
end if
!endpar
#endif

!III *********** division, SA, extrapolation, energy, tests *************

!III.1 division by denominators (divided amplitudes are stored in the same files
!      like T1n and T2n

!III.1.1 div. T1n
call divt(Work(iOff),wrksize,2,mapdt13,mapit13,mapddp1,mapidp1,mapddp2,mapidp2,rc)
call divt(Work(iOff),wrksize,2,mapdt14,mapit14,mapddp1,mapidp1,mapddp2,mapidp2,rc)

!III.1.2 div. T2n
call divt(Work(iOff),wrksize,4,mapdt21,mapit21,mapddp1,mapidp1,mapddp2,mapidp2,rc)
call divt(Work(iOff),wrksize,4,mapdt22,mapit22,mapddp1,mapidp1,mapddp2,mapidp2,rc)
call divt(Work(iOff),wrksize,4,mapdt23,mapit23,mapddp1,mapidp1,mapddp2,mapidp2,rc)

!III.2 Spin adaptation

if (keysa > 0) then
  call saamp(Work(iOff),wrksize,keysa)
end if

!III.3  extrapolation
if (yesext == 1) then

  !* write new Tn (T21,T22,T23,T13,T14) to T stack
  call diiswa1(Work(iOff),wrksize,diispointt)

  !* calc Tn = Tn - E
  call calcr(Work(iOff),wrksize,lune)

  !* write Tn=Tn-E to R stack
  call diiswa1(Work(iOff),wrksize,diispointr)

  !* do diis Tn = DIIS (Tprev) if necc.
  call diis(Work(iOff),wrksize,diispointt,diispointr,keyext)

  !* write Tn = T(DIIS) to E
  ! rewind lune
  call filemanager(2,lune,rc)
  ! T2aaaa
  call wrtmediate(Work(iOff),wrksize,lune,mapdt21,mapit21,rc)
  ! T2bbbb
  call wrtmediate(Work(iOff),wrksize,lune,mapdt22,mapit22,rc)
  ! T2abab
  call wrtmediate(Work(iOff),wrksize,lune,mapdt23,mapit23,rc)
  ! T1aa
  call wrtmediate(Work(iOff),wrksize,lune,mapdt13,mapit13,rc)
  ! T1bb
  call wrtmediate(Work(iOff),wrksize,lune,mapdt14,mapit14,rc)

end if

!III.4 save restart informations - amplitudes
if (keyrst /= 0) then
  call saverest1(Work(iOff),wrksize,lunrst)
end if

!III.5 write Tn into place of to

!III.5.1 put t1naa -> t1oaa
call map(Work(iOff),wrksize,2,1,2,0,0,mapdt13,mapit13,1,mapdt11,mapit11,posst110,posst,rc)

!III.5.2 put t1nbb -> t1obb
call map(Work(iOff),wrksize,2,1,2,0,0,mapdt14,mapit14,1,mapdt12,mapit12,posst120,posst,rc)

!III.5.3 rewind lunt2o1 and write t2naaaa there
call filemanager(2,lunt2o1,rc)
call wrtmediate(Work(iOff),wrksize,lunt2o1,mapdt21,mapit21,rc)

!III.5.4 rewind lunt2o2 and write t2nbbbb there
call filemanager(2,lunt2o2,rc)
call wrtmediate(Work(iOff),wrksize,lunt2o2,mapdt22,mapit22,rc)

!III.5.5 rewind lunt2o3 and write t2nabab there
call filemanager(2,lunt2o3,rc)
call wrtmediate(Work(iOff),wrksize,lunt2o3,mapdt23,mapit23,rc)

! test zero
call percentzero(Work(iOff),wrksize,mapdt13,pz1aa)
call percentzero(Work(iOff),wrksize,mapdt14,pz1bb)
!call percentzero(Work(iOff),wrksize,mapdt21,pz2aaaa)
!call percentzero(Work(iOff),wrksize,mapdt22,pz2bbbb)
call percentzero(Work(iOff),wrksize,mapdt23,pz2abab)

!III.6 calc energy
!FUE energy = 0.0d0
energy = Escf
scalar = 0.0d0

!3.6.1calc fai(a,i)aa . t1(a,i)aa
call multdot(Work(iOff),wrksize,2,mapdfk3,mapifk3,1,mapdt13,mapit13,1,scalar,rc)
energy = energy+scalar
E1aa = scalar
!FUE write(6,22) scalar

!III.6.2 calc fai(a,i)bb . t1(a,i)bb
call multdot(Work(iOff),wrksize,2,mapdfk4,mapifk4,1,mapdt14,mapit14,1,scalar,rc)
energy = energy+scalar
E1bb = scalar
!FUE write(6,22) scalar

!III.6.3 get <ab||ij>aaaa into V1 and calc V1(abij) . Tau(abij)aaaa
call filemanager(2,lunabij1,rc)
call getmediate(Work(iOff),wrksize,lunabij1,possv10,mapdv1,mapiv1,rc)
call mktau(Work(iOff),wrksize,mapdt21,mapit21,mapdt13,mapit13,mapdt13,mapit13,1.0d0,rc)
call multdot(Work(iOff),wrksize,4,mapdv1,mapiv1,1,mapdt21,mapit21,1,scalar,rc)
energy = energy+scalar
E2aaaa = scalar
!FUE write(6,22) scalar

!III.6.4 get <ab||ij>bbbb into V1 and calc V1(abij) . Tau(abij)bbbb
call filemanager(2,lunabij2,rc)
call getmediate(Work(iOff),wrksize,lunabij2,possv10,mapdv1,mapiv1,rc)
call mktau(Work(iOff),wrksize,mapdt22,mapit22,mapdt14,mapit14,mapdt14,mapit14,1.0d0,rc)
call multdot(Work(iOff),wrksize,4,mapdv1,mapiv1,1,mapdt22,mapit22,1,scalar,rc)
energy = energy+scalar
E2bbbb = scalar
!FUE write(6,22) scalar

!III.6.5 get <ab||ij>abab into V1 and calc V1(abij) . Tau(abij)abab
call filemanager(2,lunabij3,rc)
call getmediate(Work(iOff),wrksize,lunabij3,possv10,mapdv1,mapiv1,rc)
call mktau(Work(iOff),wrksize,mapdt23,mapit23,mapdt13,mapit13,mapdt14,mapit14,1.0d0,rc)
call multdot(Work(iOff),wrksize,4,mapdv1,mapiv1,1,mapdt23,mapit23,1,scalar,rc)
energy = energy+scalar
E2abab = scalar
!FUE write(6,22) scalar
!FUE write(6,22) energy
!22 format (2x,f20.15)

!III.7  save restart informations - energy, iteration cycle
if (keyrst /= 0) then
  call saverest2(lunrst,energy,niter,iokey,daddr(lunrst))
end if

!tmp call timing(timtotcpu,timdifcpu,timtotwc,timdifwc)
call CWTime(timtotcpun,timtotwcn)
!tmp timdifcpu = timtotcpun-timtotcpu
timdifwc = timtotwcn-timtotwc
timtotcpu = timtotcpun
timtotwc = timtotwcn
if (fullprint >= 1) then
  write(6,*) ' Adapt., Extrap. and Energy evaluation done ',myRank,timdifwc,timtotwc
end if

!III.8  test of convergence or termination

diff = energy-energyold
if (niter == 1) diff = E1aa+E1bb+E2aaaa+E2bbbb+E2abab
energyold = energy

if ((fullprint >= 0) .and. (fullprint <= 1)) then
  ! reduced printing
  write(6,'(6X,I5,5X,F17.8,2X,2(F15.8,2X))') niter,energy,E1aa+E1bb+E2aaaa+E2bbbb+E2abab,diff

else if (fullprint >= 2) then
  ! full printing
  write(6,*)
  write(6,'(6X,A,I4)') 'Iteration No        :',niter
  write(6,'(6X,A,2F24.13)') 'Total energy (diff) :',energy,diff
  write(6,'(6X,A,F24.13)') 'Correlation energy  :',E1aa+E1bb+E2aaaa+E2bbbb+E2abab
  write(6,'(6X,A,F24.13)') 'Reference energy    :',Escf
  write(6,'(6X,A,F17.8)') 'E1aa   contribution :',E1aa
  write(6,'(6X,A,F17.8)') 'E1bb   contribution :',E1bb
  write(6,'(6X,A,F17.8)') 'E2aaaa contribution :',E2aaaa
  write(6,'(6X,A,F17.8)') 'E2bbbb contribution :',E2bbbb
  write(6,'(6X,A,F17.8)') 'E2abab contribution :',E2abab
  write(6,'(6X,A)') '% of small amplitudes in  T1aa   T1bb   T2aaaa T2bbbb T2abab'
  write(6,'(31X,5(F5.1,2X))') pz1aa,pz1bb,pz2abab

end if
call Add_Info('E_CCSD',[Energy],1,8)

niter = niter+1

!--- the original code depends on print level!!!!!!!
!if ((abs(diff) <= ccconv) .and. (fullprint >= 0)) then
!  write(6,*) '     Convergence after ',niter,' Iterations'
!else if ((niter > maxiter) .and. (fullprint >= 0)) then
!  write(6,*) '     Convergence not reached'
!else if (niter <= maxiter)then
!  goto 1
!end if
if (abs(diff) <= ccconv) then
  write(6,*) '     Convergence after ',niter,' Iterations'
else if (niter > maxiter) then
  write(6,*) '     Convergence not reached'
else if (niter <= maxiter) then
  goto 1
end if
if (fullprint >= 0) then
  write(6,*)
  write(6,*)

  !* **************     Finito      **************

  write(6,'(6X,A,2F17.8)') 'Total energy (diff) :',energy,diff
  write(6,'(6X,A,F24.13)') 'Correlation energy  :',E1aa+E1bb+E2aaaa+E2bbbb+E2abab
  write(6,'(6X,A,F24.13)') 'Reference energy    :',Escf
  write(6,'(6X,A,F17.8)') 'E1aa   contribution :',E1aa
  write(6,'(6X,A,F17.8)') 'E1bb   contribution :',E1bb
  write(6,'(6X,A,F17.8)') 'E2aaaa contribution :',E2aaaa
  write(6,'(6X,A,F17.8)') 'E2bbbb contribution :',E2bbbb
  write(6,'(6X,A,F17.8)') 'E2abab contribution :',E2abab
  write(6,*)
  write(6,*)
end if
! Export a method and energy to the MOLCAS runfile
call Put_cArray('Relax Method','CCSDT   ',8)
call Store_Energies(1,[Energy],1)

!IV.0  type 5 maximal elements + euclidian norms in each type of amplitudes

!IV.0.1T1aa
call max5(Work(iOff),wrksize,2,mapdt13,mapit13,'T1aa    ')

!IV.0.2T1bb
call max5(Work(iOff),wrksize,2,mapdt14,mapit14,'T1bb    ')

!IV.0.3T2aaaa
call max5(Work(iOff),wrksize,4,mapdt21,mapit21,'T2aaaa  ')

!IV.0.4T2bbbb
call max5(Work(iOff),wrksize,4,mapdt22,mapit22,'T2bbbb  ')

!IV.0.5T2abab
call max5(Work(iOff),wrksize,4,mapdt23,mapit23,'T2abab  ')

!IV.1  close lunabij1,2,3
call filemanager(3,lunabij1,rc)
call filemanager(3,lunabij2,rc)
call filemanager(3,lunabij3,rc)

!IV.2  close lunt2o1,2,3
call filemanager(3,lunt2o1,rc)
call filemanager(3,lunt2o2,rc)
call filemanager(3,lunt2o3,rc)

!IV.3  close lunext files if extrapolation is used
if (yesext == 1) then
  call diiscf(diispointt,cycext)
  call diiscf(diispointr,cycext)
  call filemanager(3,lune,rc)
end if

!IV.4  close lunrst file if restart informations was stored
if (keyrst /= 0) then
  call filemanager(3,lunrst,rc)
end if
! global ccsd synchronization point for parallel runs
call GASync()

! Releasing the memory
call GetMem('CCSD','Free','Real',iOff,wrksize)

999 continue
if (fullprint >= 0) then
  write(6,*)
  write(6,'(6X,A)') 'Happy Landing!'
  write(6,*)
end if

ireturn = 0

return

end subroutine ccsd
