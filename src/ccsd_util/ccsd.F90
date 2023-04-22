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

use ccsd_global, only: ccconv, cycext, daddr, dp1, dp2, Escf, fk1, fk2, fk3, fk4, fk5, fk6, fullprint, idaaaa, idaabb, idab, &
                       idabba, idbaab, idbbaa, idbbbb, ideffab, ididle, idfin, idtmab, iokey, keyrst, keysa, maxiter, maxspace, n, &
                       noccsd, nprocab, p, t11, t12, t13, t14, t21, t22, t23, v1, v2, v3, v4, yesext
use Para_Info, only: MyRank, nProcs
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
logical(kind=iwp), intent(inout) :: run_triples
integer(kind=iwp) :: diispointr(4), diispointt(4), i, idum(1), infree, inv4, istatus, keyexc, keyext, lenn, lenv, lunabij1, &
                     lunabij2, lunabij3, lune, lunrst, lunt2o1, lunt2o2, lunt2o3, lunw3aaaa, lunw3aabb, lunw3abba, lunw3baab, &
                     lunw3bbaa, lunw3bbbb, nabstack, nfree, niter, posabstack, post, rc, wrksize
real(kind=wp) :: diff, dum(1), E1aa, E1bb, E2aaaa, E2abab, E2bbbb, energy, energyold, pz1aa, pz1bb, pz2abab, scalar, timdifwc, &
                 timtotcpu, timtotcpun, timtotit, timtotwc, timtotwcn
real(kind=wp), allocatable :: wrk(:)
integer(kind=iwp), external :: iPrintLevel

call CWTime(timtotcpu,timtotwc)
timtotit = timtotwc

energy = 9999999.0_wp
energyold = 99999.0_wp
fullprint = 0
if (iPrintLevel(-1) <= 0) fullprint = -1

!I **************  start   section **************

!I.1.1 read input data form INPDAT (reorg) and input file

call reainput()

!I.1.2 write head to output file

if (fullprint >= 0) then
  call wrhead()
  if (noccsd == 1) then
    write(u6,'(6X,A)') 'NOSD key is turned on'
    write(u6,'(6X,A)') 'CCSD step is skipped'
    call happy()
    return
  end if
  write(u6,*) ' nProcs, myRank',nProcs,myRank
end if

!I.1.3 check if there is something to do in CCSD
call ccsd_exc(keyexc)
if (keyexc == 0) then
  write(u6,*) ' No determinants in CCSD expansion '
  write(u6,*) ' Zero correlation energy           '
  write(u6,*) ' CCSD step excluded                '
  if (run_triples) then
    write(u6,*) ' Triples step excluded'
    run_triples = .false.
  end if
  call happy()
  return
else if (keyexc == 1) then
  write(u6,*) ' Only monoexcitations in CC expansion '
  write(u6,*) ' CCSD step excluded                '
  if (run_triples) then
    write(u6,*) ' Triples step excluded'
    run_triples = .false.
  end if
  call happy()
  return
end if

!I.1.4 def parameters for filemgr

call mkfilemgrcom()

!I.2.1 calc. work space requirements to fix and help mediates

call initfiles(wrksize,lenv,lenn)
if (fullprint >= 0) write(u6,*) ' Basic Work space requirements    :',wrksize

!I.2.3 allocate work space

!I.2.3.1 check free space
call mma_maxDBLE(maxspace)
maxspace = maxspace-8
if (fullprint >= 0) write(u6,*) ' Max Size',maxspace

if (maxspace < wrksize) then
  write(u6,*) ' Allocation of work space failed'
  write(u6,*) ' Increase the value of the variable MOLCAS_MEM'
  call Abend()
end if

!I.2.3.2 fix, where AB stack will be located
nfree = maxspace-wrksize
inv4 = int(lenv/lenn)
infree = int(nfree/lenn)
if (infree > int(inv4/2)) then
  ! there is enough free room for AB_stack
  posabstack = wrksize+1
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
  posabstack = v4%pos0+nabstack*lenn
end if

if (fullprint >= 2) write(u6,*) ' Dimension of AB stack            :',nabstack
if (fullprint >= 0) write(u6,*) ' Final Work space requirements    :',wrksize

!I.2.3.3 check free space
call mma_allocate(wrk,wrksize,label='CCSD')

!I.2.4 set wrk = 0
wrk(:) = Zero
if (fullprint >= 0) write(u6,*) ' Allocation of work space   : Done'

!I.2.5 read static integrals from INTSTA (reorg) file

call reaintsta(wrk,wrksize)

!I.2.6 divide fok to faa,fai,fii and dp

call divfok(wrk,wrksize,n,p,fk1,fk2,fk3,fk4,fk5,fk6,dp1,dp2,rc)

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
  call getmediate(wrk,wrksize,lunrst,t11,rc)
  ! get T1bb
  call getmediate(wrk,wrksize,lunrst,t12,rc)

  ! get and write T2aaaa
  call getmediate(wrk,wrksize,lunrst,v4,rc)
  call wrtmediate(wrk,wrksize,lunt2o1,v4,rc)
  ! get and write T2bbbb
  call getmediate(wrk,wrksize,lunrst,v4,rc)
  call wrtmediate(wrk,wrksize,lunt2o2,v4,rc)
  ! get and write T2abab
  call getmediate(wrk,wrksize,lunrst,v4,rc)
  call wrtmediate(wrk,wrksize,lunt2o3,v4,rc)

  if (iokey == 1) then
    ! Fortran IO
    read(lunrst,iostat=istatus) energyold,niter
    if (istatus < 0) then
      ! case when energy and niter is not in save file
      !FUE energyold = Zero
      energyold = Escf
      niter = 0
      write(u6,*) ' ENERGY AND NIT WAS NOT IN SAVE FILE, CHANGED TO 0'
    end if

  else
    ! MOLCAS IO
    call ddafile(lunrst,2,dum,1,daddr(lunrst))
    energyold = dum(1)
    call idafile(lunrst,2,idum,1,daddr(lunrst))
    niter = idum(1)
  end if

  niter = niter+1

else
  !I.2.9.3 T2naaaa,T2nbbbb and T2nabab are 0, we can write them as T2old
  call wrtmediate(wrk,wrksize,lunt2o1,t21,rc)
  call wrtmediate(wrk,wrksize,lunt2o2,t22,rc)
  call wrtmediate(wrk,wrksize,lunt2o3,t23,rc)
end if

!I.2.10 write V1 (<ab||ij>aaaa) to lunabij1
! V2 (<ab||ij>bbbb) to lunabij2
! V3 (<ab||ij>abab) to lunabij3

call wrtmediate(wrk,wrksize,lunabij1,v1,rc)
call wrtmediate(wrk,wrksize,lunabij2,v2,rc)
call wrtmediate(wrk,wrksize,lunabij3,v3,rc)

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
  call wrtmediate(wrk,wrksize,lune,t21,rc)
  ! T2bbbb
  call wrtmediate(wrk,wrksize,lune,t22,rc)
  ! T2abab
  call wrtmediate(wrk,wrksize,lune,t23,rc)
  ! T1aa
  call wrtmediate(wrk,wrksize,lune,t11,rc)
  ! T1bb
  call wrtmediate(wrk,wrksize,lune,t12,rc)

  keyext = 0
end if

!* write head for short + medium printing
if ((fullprint >= 0) .and. (fullprint < 2)) then
  write(u6,*)
  write(u6,'(6X,A)') 'Iteration       Total enegy      Corr. energy      Difference'
end if

!II **************  Iteration cycle   **************

if (keyrst /= 2) niter = 1

call distnodes()

do
  if ((fullprint >= 1) .and. (myRank == 0) .and. (nProcs > 1)) then
    write(u6,*) 'Npocab:',nprocab
    do i=1,nprocab
      write(u6,*) 'Proc  :',i,idab(i),ideffab(i)
    end do
    write(u6,*) 'IDaaaa:',idaaaa
    write(u6,*) 'IDaabb:',idaabb
    write(u6,*) 'IDabba:',idabba
    write(u6,*) 'IDbbbb:',idbbbb
    write(u6,*) 'IDbbaa:',idbbaa
    write(u6,*) 'IDbaab:',idbaab
    write(u6,*) 'IDfin :',idfin
  end if
  !par
  ! set ididle,idtmab =0
  ididle(1:nProcs) = Zero
  idtmab(1:nProcs) = Zero
  !end par

  call init(wrk,wrksize,lunabij1,lunabij2,lunabij3)
  call CWTime(timtotcpun,timtotwcn)
  !tmp timdifcpu = timtotcpun-timtotcpu
  timdifwc = timtotwcn-timtotwc
  timtotcpu = timtotcpun
  timtotwc = timtotwcn
  if (fullprint >= 1) then
    write(u6,*) ' Initial section done',myRank,timdifwc,timtotwcn-timtotit
    timtotit = timtotwcn
  end if

  call sumoverab(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3,nabstack,posabstack)
  call CWTime(timtotcpun,timtotwcn)
  !tmp timdifcpu = timtotcpun-timtotcpu
  timdifwc = timtotwcn-timtotwc
  timtotcpu = timtotcpun
  timtotwc = timtotwcn
  if (fullprint >= 1) write(u6,*) ' Summation over ab done',myRank,timdifwc

  !par store time spend in sumoverab process for this node
  idtmab(myRank+1) = timdifwc

  call sumovera(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3,lunw3aaaa,lunw3baab,lunw3bbaa,lunw3bbbb,lunw3abba,lunw3aabb)
  call CWTime(timtotcpun,timtotwcn)
  !tmp timdifcpu = timtotcpun-timtotcpu
  timdifwc = timtotwcn-timtotwc
  timtotcpu = timtotcpun
  timtotwc = timtotwcn
  if (fullprint >= 1) write(u6,*) ' Summation over a done',myRank,timdifwc

  call intermezzo(wrk,wrksize,lunw3aaaa,lunw3bbbb,lunw3abba,lunw3baab,lunw3aabb,lunw3bbaa,lunt2o1,lunt2o2,lunt2o3,lunabij1, &
                  lunabij2,lunabij3)
  call CWTime(timtotcpun,timtotwcn)
  !tmp timdifcpu = timtotcpun-timtotcpu
  timdifwc = timtotwcn-timtotwc
  timtotcpu = timtotcpun
  timtotwc = timtotwcn
  if (fullprint >= 1) write(u6,*) ' Internal section done',myRank,timdifwc

  call finale(wrk,wrksize,lunabij1,lunabij2,lunabij3,lunt2o1,lunt2o2,lunt2o3)
  call CWTime(timtotcpun,timtotwcn)
  !tmp timdifcpu = timtotcpun-timtotcpu
  timdifwc = timtotwcn-timtotwc
  timtotcpu = timtotcpun
  timtotwc = timtotwcn
  if (fullprint >= 1) write(u6,*) ' Final section done',myRank,timdifwc

# ifdef _MOLCAS_MPP_
  !par
  call joinamplitudes(wrk,wrksize)
  call CWTime(timtotcpun,timtotwcn)
  !tmp timdifcpu = timtotcpun-timtotcpu
  timdifwc = timtotwcn-timtotwc
  timtotcpu = timtotcpun
  timtotwc = timtotwcn
  if (fullprint >= 1) write(u6,*) ' Amplitudes allreduced',myRank,timdifwc

  ! store differential wc time at this point to calculate idle time
  ! and redefine ideffab for next step
  ididle(myRank+1) = timdifwc
  if (niter > 0) then
    call redef()
    if ((myRank == 0) .and. (fullprint >= 1)) then
      write(u6,997) ' ABTim',(idtmab(i),i=1,nProcs)
      write(u6,997) ' IDtim',(ididle(i),i=1,nProcs)
      write(u6,997) ' EFF  ',(ideffab(i),i=1,nprocab)
    end if
  end if
  !endpar
# endif

  !III *********** division, SA, extrapolation, energy, tests *************

  !III.1 division by denominators (divided amplitudes are stored in the same files
  !      like T1n and T2n

  !III.1.1 div. T1n
  call divt(wrk,wrksize,2,t13,dp1,dp2,rc)
  call divt(wrk,wrksize,2,t14,dp1,dp2,rc)

  !III.1.2 div. T2n
  call divt(wrk,wrksize,4,t21,dp1,dp2,rc)
  call divt(wrk,wrksize,4,t22,dp1,dp2,rc)
  call divt(wrk,wrksize,4,t23,dp1,dp2,rc)

  !III.2 Spin adaptation

  if (keysa > 0) call saamp(wrk,wrksize,keysa)

  !III.3  extrapolation
  if (yesext == 1) then

    !* write new Tn (T21,T22,T23,T13,T14) to T stack
    call diiswa1(wrk,wrksize,diispointt)

    !* calc Tn = Tn - E
    call calcr(wrk,wrksize,lune)

    !* write Tn=Tn-E to R stack
    call diiswa1(wrk,wrksize,diispointr)

    !* do diis Tn = DIIS (Tprev) if necc.
    call diis(wrk,wrksize,diispointt,diispointr,keyext)

    !* write Tn = T(DIIS) to E
    ! rewind lune
    call filemanager(2,lune,rc)
    ! T2aaaa
    call wrtmediate(wrk,wrksize,lune,t21,rc)
    ! T2bbbb
    call wrtmediate(wrk,wrksize,lune,t22,rc)
    ! T2abab
    call wrtmediate(wrk,wrksize,lune,t23,rc)
    ! T1aa
    call wrtmediate(wrk,wrksize,lune,t13,rc)
    ! T1bb
    call wrtmediate(wrk,wrksize,lune,t14,rc)

  end if

  !III.4 save restart informations - amplitudes
  if (keyrst /= 0) call saverest1(wrk,wrksize,lunrst)

  !III.5 write Tn into place of to

  !III.5.1 put t1naa -> t1oaa
  call map(wrk,wrksize,2,1,2,0,0,t13,1,t11,post,rc)

  !III.5.2 put t1nbb -> t1obb
  call map(wrk,wrksize,2,1,2,0,0,t14,1,t12,post,rc)

  !III.5.3 rewind lunt2o1 and write t2naaaa there
  call filemanager(2,lunt2o1,rc)
  call wrtmediate(wrk,wrksize,lunt2o1,t21,rc)

  !III.5.4 rewind lunt2o2 and write t2nbbbb there
  call filemanager(2,lunt2o2,rc)
  call wrtmediate(wrk,wrksize,lunt2o2,t22,rc)

  !III.5.5 rewind lunt2o3 and write t2nabab there
  call filemanager(2,lunt2o3,rc)
  call wrtmediate(wrk,wrksize,lunt2o3,t23,rc)

  ! test zero
  call percentzero(wrk,wrksize,t13,pz1aa)
  call percentzero(wrk,wrksize,t14,pz1bb)
  !call percentzero(wrk,wrksize,t21,pz2aaaa)
  !call percentzero(wrk,wrksize,t22,pz2bbbb)
  call percentzero(wrk,wrksize,t23,pz2abab)

  !III.6 calc energy
  !FUE energy = Zero
  energy = Escf
  scalar = Zero

  !3.6.1calc fai(a,i)aa . t1(a,i)aa
  call multdot(wrk,wrksize,2,fk3,1,t13,1,scalar,rc)
  energy = energy+scalar
  E1aa = scalar
  !FUE write(u6,22) scalar

  !III.6.2 calc fai(a,i)bb . t1(a,i)bb
  call multdot(wrk,wrksize,2,fk4,1,t14,1,scalar,rc)
  energy = energy+scalar
  E1bb = scalar
  !FUE write(u6,22) scalar

  !III.6.3 get <ab||ij>aaaa into V1 and calc V1(abij) . Tau(abij)aaaa
  call filemanager(2,lunabij1,rc)
  call getmediate(wrk,wrksize,lunabij1,v1,rc)
  call mktau(wrk,wrksize,t21,t13,t13,One,rc)
  call multdot(wrk,wrksize,4,v1,1,t21,1,scalar,rc)
  energy = energy+scalar
  E2aaaa = scalar
  !FUE write(u6,22) scalar

  !III.6.4 get <ab||ij>bbbb into V1 and calc V1(abij) . Tau(abij)bbbb
  call filemanager(2,lunabij2,rc)
  call getmediate(wrk,wrksize,lunabij2,v1,rc)
  call mktau(wrk,wrksize,t22,t14,t14,One,rc)
  call multdot(wrk,wrksize,4,v1,1,t22,1,scalar,rc)
  energy = energy+scalar
  E2bbbb = scalar
  !FUE write(u6,22) scalar

  !III.6.5 get <ab||ij>abab into V1 and calc V1(abij) . Tau(abij)abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v1,rc)
  call mktau(wrk,wrksize,t23,t13,t14,One,rc)
  call multdot(wrk,wrksize,4,v1,1,t23,1,scalar,rc)
  energy = energy+scalar
  E2abab = scalar
  !FUE write(u6,22) scalar
  !FUE write(u6,22) energy

  !III.7  save restart informations - energy, iteration cycle
  if (keyrst /= 0) call saverest2(lunrst,energy,niter,iokey,daddr(lunrst))

  !tmp call timing(timtotcpu,timdifcpu,timtotwc,timdifwc)
  call CWTime(timtotcpun,timtotwcn)
  !tmp timdifcpu = timtotcpun-timtotcpu
  timdifwc = timtotwcn-timtotwc
  timtotcpu = timtotcpun
  timtotwc = timtotwcn
  if (fullprint >= 1) write(u6,*) ' Adapt., Extrap. and Energy evaluation done ',myRank,timdifwc,timtotwc

  !III.8  test of convergence or termination

  diff = energy-energyold
  if (niter == 1) diff = E1aa+E1bb+E2aaaa+E2bbbb+E2abab
  energyold = energy

  if ((fullprint >= 0) .and. (fullprint <= 1)) then
    ! reduced printing
    write(u6,'(6X,I5,5X,F17.8,2X,2(F15.8,2X))') niter,energy,E1aa+E1bb+E2aaaa+E2bbbb+E2abab,diff

  else if (fullprint >= 2) then
    ! full printing
    write(u6,*)
    write(u6,'(6X,A,I4)') 'Iteration No        :',niter
    write(u6,'(6X,A,2F24.13)') 'Total energy (diff) :',energy,diff
    write(u6,'(6X,A,F24.13)') 'Correlation energy  :',E1aa+E1bb+E2aaaa+E2bbbb+E2abab
    write(u6,'(6X,A,F24.13)') 'Reference energy    :',Escf
    write(u6,'(6X,A,F17.8)') 'E1aa   contribution :',E1aa
    write(u6,'(6X,A,F17.8)') 'E1bb   contribution :',E1bb
    write(u6,'(6X,A,F17.8)') 'E2aaaa contribution :',E2aaaa
    write(u6,'(6X,A,F17.8)') 'E2bbbb contribution :',E2bbbb
    write(u6,'(6X,A,F17.8)') 'E2abab contribution :',E2abab
    write(u6,'(6X,A)') '% of small amplitudes in  T1aa   T1bb   T2aaaa T2bbbb T2abab'
    write(u6,'(31X,5(F5.1,2X))') pz1aa,pz1bb,pz2abab

  end if
  call Add_Info('E_CCSD',[Energy],1,8)

  niter = niter+1

  if (abs(diff) <= ccconv) then
    write(u6,*) '     Convergence after ',niter,' Iterations'
    exit
  else if (niter > maxiter) then
    write(u6,*) '     Convergence not reached'
    exit
  end if
end do
if (fullprint >= 0) then
  write(u6,*)
  write(u6,*)

  !* **************     Finito      **************

  write(u6,'(6X,A,2F17.8)') 'Total energy (diff) :',energy,diff
  write(u6,'(6X,A,F24.13)') 'Correlation energy  :',E1aa+E1bb+E2aaaa+E2bbbb+E2abab
  write(u6,'(6X,A,F24.13)') 'Reference energy    :',Escf
  write(u6,'(6X,A,F17.8)') 'E1aa   contribution :',E1aa
  write(u6,'(6X,A,F17.8)') 'E1bb   contribution :',E1bb
  write(u6,'(6X,A,F17.8)') 'E2aaaa contribution :',E2aaaa
  write(u6,'(6X,A,F17.8)') 'E2bbbb contribution :',E2bbbb
  write(u6,'(6X,A,F17.8)') 'E2abab contribution :',E2abab
  write(u6,*)
  write(u6,*)
end if
! Export a method and energy to the MOLCAS runfile
call Put_cArray('Relax Method','CCSDT   ',8)
call Store_Energies(1,[Energy],1)

!IV.0 type 5 maximal elements + euclidian norms in each type of amplitudes

!IV.0.1 T1aa
call max5(wrk,wrksize,2,t13,'T1aa    ')

!IV.0.2 T1bb
call max5(wrk,wrksize,2,t14,'T1bb    ')

!IV.0.3 T2aaaa
call max5(wrk,wrksize,4,t21,'T2aaaa  ')

!IV.0.4 T2bbbb
call max5(wrk,wrksize,4,t22,'T2bbbb  ')

!IV.0.5 T2abab
call max5(wrk,wrksize,4,t23,'T2abab  ')

!IV.1 close lunabij1,2,3
call filemanager(3,lunabij1,rc)
call filemanager(3,lunabij2,rc)
call filemanager(3,lunabij3,rc)

!IV.2 close lunt2o1,2,3
call filemanager(3,lunt2o1,rc)
call filemanager(3,lunt2o2,rc)
call filemanager(3,lunt2o3,rc)

!IV.3 close lunext files if extrapolation is used
if (yesext == 1) then
  call diiscf(diispointt,cycext)
  call diiscf(diispointr,cycext)
  call filemanager(3,lune,rc)
end if

!IV.4  close lunrst file if restart informations was stored
if (keyrst /= 0) call filemanager(3,lunrst,rc)
! global ccsd synchronization point for parallel runs
call GASync()

! Releasing the memory
call mma_deallocate(wrk)

call happy()

return

!FUE 22 format(2x,f20.15)
# ifdef _MOLCAS_MPP_
997 format(a6,2x,16(f7.2,1x))
#endif

contains

subroutine happy()

  if (fullprint >= 0) then
    write(u6,*)
    write(u6,'(6X,A)') 'Happy Landing!'
    write(u6,*)
  end if

  ireturn = 0

end subroutine happy

end subroutine ccsd
