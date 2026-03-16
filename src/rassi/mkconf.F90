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

subroutine MKCONF(ICNFTAB)

use definitions, only: iwp, u6
use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: nSym => nIrrep, MUL

implicit none
integer(kind=iwp), intent(inout) :: ICNFTAB(*)
integer(kind=iwp) NEL, MINOP, MAXOP, LSYM
integer(kind=iwp) MXPRT
parameter(MXPRT=150)
integer(kind=iwp) LIMPOP(2,MXPRT), LIMOP(2,MXPRT), LIMCL(2,MXPRT)
integer(kind=iwp) IOPDST(MXPRT), ICLDST(MXPRT), IOC(MXPRT), ICNF(MXPRT)
integer(kind=iwp) LIM1, LIM2, LIM1SUM, LIM2SUM, IGAS, NOR, IERR
integer(kind=iwp) MNOP, MXOP, MNCL, MXCL, NOPN, NCLS
integer(kind=iwp) INIT1, INIT2, M, MORE, NOP1, NCL2
integer(kind=iwp) ISYM, NORB, ICONF
integer(kind=iwp) IFORM, IR, ITYPE, IW, KCNFSTA, KGASLIM
integer(kind=iwp) KGASORB, KINFO, KPOS, LCLS, LENCNF, LOPN, NCNF
integer(kind=iwp) NCNFSYM(8), NGAS, NOCC, NTAB
integer(kind=iwp) I, J, K, IOFF, IO, IORB, N, NCL, NOP
integer(kind=iwp), allocatable :: ISM(:)
intrinsic MIN, MAX
integer(kind=iwp) :: IPOW4(0:15) = [1,4,16,64,256,1024,4096,16384,65536,262144,1048576,4194304,16777216,67108864,268435456, &
                                    1073741824]
integer(kind=iwp) :: IPOW256(0:3) = [1,256,65536,16777216]

ITYPE = ICNFTAB(2)
if (ITYPE /= 37) then
  write(u6,*) 'MKCONF error: This is not a configuration table!'
  call ABEND()
end if
! Unbutton the CNF table.
NTAB = ICNFTAB(1)
NEL = ICNFTAB(3)
NORB = ICNFTAB(4)
MINOP = ICNFTAB(5)
MAXOP = ICNFTAB(6)
NSYM = ICNFTAB(7)
LSYM = ICNFTAB(8)
NGAS = ICNFTAB(9)
IFORM = ICNFTAB(10)
KGASORB = 11
KGASLIM = KGASORB+(NSYM+1)*(NGAS+1)
! Check and refine the GAS limits:
if ((NGAS <= 0) .or. (NGAS > MXPRT)) then
  write(u6,*) ' MKCONF ERROR: Nr of GAS partitions is out of'
  write(u6,*) ' bounds. NGAS must be > 0 and < MXPRT=',MXPRT
  write(u6,*) ' Input argument NGAS is ',NGAS
  call ABEND()
end if

LIM1SUM = 0
LIM2SUM = 0
NORB = 0
do IGAS=1,NGAS
  NOR = ICNFTAB(KGASORB+(NSYM+1)*IGAS)
  if (NOR < 0) IERR = 1
  LIM1 = max(0,ICNFTAB(KGASLIM+2*(IGAS-1)))
  LIM2 = min(NEL,ICNFTAB(KGASLIM+1+2*(IGAS-1)),2*NOR)
  LIM1SUM = LIM1SUM+LIM1
  LIM2SUM = LIM2SUM+LIM2
  LIMPOP(1,IGAS) = LIM1
  LIMPOP(2,IGAS) = LIM2
  NORB = NORB+NOR
end do

IERR = 0
do IGAS=1,NGAS
  LIM1 = max(LIMPOP(1,IGAS),NEL-(LIM2SUM-LIMPOP(2,IGAS)))
  LIM2 = min(LIMPOP(2,IGAS),NEL-(LIM1SUM-LIMPOP(1,IGAS)))
  LIMPOP(1,IGAS) = LIM1
  LIMPOP(2,IGAS) = LIM2
  if (LIM1 > LIM2) IERR = 1
end do

if (IERR > 0) then
  write(u6,*) ' MKCONF ERROR: The input GAS restrictions are'
  write(u6,*) ' impossible to meet. No configurations are'
  write(u6,*) ' generated. The program stops here.'
  write(u6,'(1X,A,I2)') ' Number of GAS partitions:',NGAS
  write(u6,'(1X,A,50I3)') ' Partition:',(IGAS,IGAS=1,NGAS)
  write(u6,'(1X,A,50I3)') ' NGASORB:  ',(ICNFTAB(KGASORB+(NSYM+1)*IGAS),IGAS=1,NGAS)
  write(u6,'(1X,A,50I3)') 'NGASLIM(1):',(ICNFTAB(KGASLIM+2*(IGAS-1)),IGAS=1,NGAS)
  write(u6,'(1X,A,50I3)') 'NGASLIM(2):',(ICNFTAB(KGASLIM+1+2*(IGAS-1)),IGAS=1,NGAS)
  write(u6,'(1X,A,I2)') ' Number of electrons:',NEL
  call ABEND()
end if

if ((NEL < 0) .or. (NEL > 2*NORB)) then
  write(u6,*) ' MKCONF ERROR: Nr of electrons is out of bounds.'
  write(u6,*) ' NEL must be > 0 and < 2*NORB=',2*NORB
  write(u6,*) ' Input argument NEL is ',NEL
  call ABEND()
end if

! Array for orbital symmetry:
call mma_allocate(ISM,NORB,Label='ISM')
! Initialize table with orbital symmetry.
IORB = 0
do IGAS=1,NGAS
  do ISYM=1,NSYM
    N = ICNFTAB(KGASORB+ISYM+(NSYM+1)*IGAS)
    do K=1,N
      IORB = IORB+1
      ISM(IORB) = ISYM
    end do
  end do
end do

! INFO table inside ICNFTAB:
KINFO = KGASLIM+2*NGAS
! Note: Nr of conf, their position and length can now be accessed as:
!      NCNF   =ICNFTAB(KINFO+0+3*(ISYM-1+NSYM*(NOPN-MINOP)))
!      KCNFSTA=ICNFTAB(KINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))
!      LENCNF =ICNFTAB(KINFO+2+3*(ISYM-1+NSYM*(NOPN-MINOP)))
! Make a list of possible number of open shells in each partition:
do IGAS=1,NGAS
  NOR = ICNFTAB(KGASORB+(NSYM+1)*IGAS)
  MNOP = 1
  MXOP = 0
  do I=LIMPOP(1,IGAS),LIMPOP(2,IGAS)
    MNOP = min(MNOP,mod(I,2))
    MXOP = max(MXOP,min(I,2*NOR-I))
  end do
  LIMOP(1,IGAS) = MNOP
  LIMOP(2,IGAS) = MXOP
end do

! Counter of configurations:
ICONF = 0
! Loop over the requested range of open shells:
outer: do NOPN=MINOP,MAXOP
  NCLS = (NEL-NOPN)/2
  if (NCLS < 0) cycle OUTER
  if (2*NCLS+NOPN /= NEL) cycle OUTER
  ! Size of each entry in the configuration table:
  NOCC = NCLS+NOPN
  LENCNF = NOCC
  if (IFORM == 2) LENCNF = NORB
  if (IFORM == 3) LENCNF = (NOCC+3)/4
  if (IFORM == 4) LENCNF = (NORB+14)/15
  ! Counter of configurations/symmetry for this nr of open shells:
  do ISYM=1,NSYM
    NCNFSYM(ISYM) = 0
  end do
  ! Loop over all ways of distributing NOPN open shells among
  ! the partitions. First make a start distribution:
  INIT1 = NGAS
  NOP1 = NOPN

10 continue
  ! Create the lexically lowest distribution with NOP1 open
  ! shells among the INIT1 lowest partitions, and
  ! increment the next higher partition (if any).
  ! Let M=Max tot nr of open shells in lower partitions.
  if (INIT1 < NGAS) IOPDST(INIT1+1) = IOPDST(INIT1+1)+1
  M = NOP1
  do IGAS=INIT1,1,-1
    M = M-LIMOP(1,IGAS)
  end do
  if (M < 0) cycle OUTER

  ! But actually, we start with zero. So M open shells must be
  ! distributed in excess of the allowed minimum, among the INIT1
  ! partitions.
  do IGAS=1,INIT1
    MORE = min(LIMOP(2,IGAS)-LIMOP(1,IGAS),M)
    IOPDST(IGAS) = LIMOP(1,IGAS)+MORE
    M = M-MORE
  end do
  if (M > 0) cycle OUTER
  ! At this point of the code, all possible distributions of
  ! open shells will be generated. Use them.
  ! First, use it to generate a table of limits for the distribution
  ! of closed shells:
  do IGAS=1,NGAS
    NOP = IOPDST(IGAS)
    NOR = ICNFTAB(KGASORB+(NSYM+1)*IGAS)-NOP
    LIM1 = max(0,(LIMPOP(1,IGAS)-NOP)/2)
    LIM2 = min(NOR,(LIMPOP(2,IGAS)-NOP)/2)
    if (LIM1 > LIM2) goto 110
    MNCL = LIM2
    MXCL = LIM1
    do I=LIM1,LIM2
      N = 2*I+NOP
      if ((N >= LIMPOP(1,IGAS)) .and. (N <= LIMPOP(2,IGAS))) then
        MNCL = min(MNCL,I)
        MXCL = max(MXCL,I)
      end if
    end do
    if (MNCL > MXCL) goto 110
    LIMCL(1,IGAS) = MNCL
    LIMCL(2,IGAS) = MXCL
  end do

  ! Loop over all possible ways of distributing NCLS closed shells
  ! among  the partitions, subject to restrictions.
  ! In order to create the start distribution:
  INIT2 = NGAS
  NCL2 = (NEL-NOPN)/2

20 continue
  ! Create the lexically lowest distribution with NCL2 closed shells
  ! among the INIT2 lowest partitions, and increment the next higher
  ! partition (if any).
  if (INIT2 < NGAS) ICLDST(INIT2+1) = ICLDST(INIT2+1)+1
  M = NCL2
  do IGAS=INIT2,1,-1
    M = M-LIMCL(1,IGAS)
  end do
  if (M < 0) goto 110
  do IGAS=1,INIT2
    MORE = min(LIMCL(2,IGAS)-LIMCL(1,IGAS),M)
    ICLDST(IGAS) = LIMCL(1,IGAS)+MORE
    M = M-MORE
  end do
  if (M > 0) goto 110
  ! Here follows code to use this population distribution.
  ! Initialize the configuration subarrays of partitions nr
  ! 1..NGAS, within this population distribution:
  IORB = 0
  do IGAS=1,NGAS
    NCL = ICLDST(IGAS)
    do I=1,NCL
      IORB = IORB+1
      IOC(IORB) = 2
    end do
    NOP = IOPDST(IGAS)
    do I=1,NOP
      IORB = IORB+1
      IOC(IORB) = 1
    end do
    NOR = ICNFTAB(KGASORB+(NSYM+1)*IGAS)
    do I=1,NOR-NCL-NOP
      IORB = IORB+1
      IOC(IORB) = 0
    end do
  end do
30 continue
  ! Here finally we will get all possible configurations, restricted
  ! by the population arrays. Screening by combined symmetry:
  ISYM = 1
  do IO=1,NORB
    if (IOC(IO) == 1) ISYM = MUL(ISM(IO),ISYM)
  end do
  ! Skip if wrong symmetry:
  if (.not. ((LSYM > 0) .and. (ISYM /= LSYM))) then
    ICONF = ICONF+1
    ! Put this configuration into the ICNFTAB table.
    ! First, determine where it should go:
    N = NCNFSYM(ISYM)
    NCNFSYM(ISYM) = N+1
    KCNFSTA = ICNFTAB(KINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))
    KPOS = KCNFSTA+N*LENCNF
    if (KPOS+LENCNF-1 > NTAB) then
      write(u6,*) ' MKCONF error: Table overflow.'
      write(u6,*) 'KCNFSTA:',KCNFSTA
      write(u6,*) ' LENCNF:',LENCNF
      write(u6,*) '   KPOS:',KPOS
      write(u6,*) '   NTAB:',NTAB
      call ABEND()
    end if
    ! Put together configuration array in standard format:
    LCLS = 1
    LOPN = NCLS+1
    do IO=1,NORB
      N = IOC(IO)
      if (N == 1) then
        ICNF(LOPN) = IO
        LOPN = LOPN+1
      else if (N == 2) then
        ICNF(LCLS) = IO
        LCLS = LCLS+1
      end if
    end do
    ! Add this configuration to the configuration table:
    if (IFORM == 1) then
      do I=1,NOCC
        ICNFTAB(KPOS-1+I) = ICNF(I)
      end do
    else if (IFORM == 3) then
      do I=1,NOCC
        IW = (3+I)/4
        IR = (3+I)-4*IW
        if (IR == 0) then
          ICNFTAB(KPOS-1+IW) = ICNF(I)
        else
          ICNFTAB(KPOS-1+IW) = ICNFTAB(KPOS-1+IW)+IPOW256(IR)*ICNF(I)
        end if
      end do
    else
      if (IFORM == 2) then
        do I=1,NORB
          ICNFTAB(KPOS-1+I) = IOC(I)
        end do
      else if (IFORM == 4) then
        do I=1,NORB
          IW = (14+I)/15
          IR = (14+I)-15*IW
          if (IR == 0) then
            ICNFTAB(KPOS-1+IW) = IOC(I)
          else
            ICNFTAB(KPOS-1+IW) = ICNFTAB(KPOS-1+IW)+IPOW4(IR)*IOC(I)

          end if
        end do
      end if
    end if

  end if
  ! Get next configuration.
  IOFF = 0
  do IGAS=1,NGAS
    NOR = ICNFTAB(KGASORB+(NSYM+1)*IGAS)
    ! Try to find next permutation within this partition:
    do K=2,NOR
      if (IOC(IOFF+K-1) > IOC(IOFF+K)) then
        do I=1,(K-1)/2
          J = IOC(IOFF+I)
          IOC(IOFF+I) = IOC(IOFF+K-I)
          IOC(IOFF+K-I) = J
        end do
        do I=K-1,1,-1
          if (IOC(IOFF+I) > IOC(IOFF+K)) then
            J = IOC(IOFF+I)
            IOC(IOFF+I) = IOC(IOFF+K)
            IOC(IOFF+K) = J
            ! OK, the next permutation has been obtained.
            goto 30
          end if
        end do
      end if
    end do
    ! Not possible. Reset permutation in this partition, and
    ! then try the next one:
    N = ICLDST(IGAS)
    M = IOPDST(IGAS)
    do IO=1,N
      IOC(IOFF+IO) = 2
    end do
    do IO=N+1,N+M
      IOC(IOFF+IO) = 1
    end do
    do IO=N+M+1,NOR
      IOC(IOFF+IO) = 0
    end do
    IOFF = IOFF+NOR
  end do
  ! All failed. No more configurations with this distribution of
  ! closed and open shells.
  ! Next ICLDST distribution. First find the first increasable index:
  M = 0
  NCL2 = -1
  do IGAS=1,NGAS
    INIT2 = IGAS-1
    if ((M > 0) .and. (ICLDST(IGAS) < LIMCL(2,IGAS))) goto 20
    M = M+ICLDST(IGAS)-LIMCL(1,IGAS)
    NCL2 = NCL2+ICLDST(IGAS)
  end do

  ! No more ICLDST distribution is possible.
110 continue

  ! Next IOPDST distribution. First find the first increasable index:
  ! That is the first partition with less than LIMOP(2,IGAS) open
  ! shells, above partitions with nonzero excess number M.
  M = 0
  NOP1 = -1
  do IGAS=1,NGAS
    INIT1 = IGAS-1
    if ((M > 0) .and. (IOPDST(IGAS) < LIMOP(2,IGAS))) goto 10
    M = M+IOPDST(IGAS)-LIMOP(1,IGAS)
    NOP1 = NOP1+IOPDST(IGAS)
  end do
  ! Temporary check: Has everything worked perfectly??
  IERR = 0
  do ISYM=1,NSYM
    N = ICNFTAB(KINFO+3*(ISYM-1+NSYM*(NOPN-MINOP)))
    if (NCNFSYM(ISYM) /= N) IERR = 1
  end do
  if (IERR /= 0) call ErrorTrap()
  ! No more IOPDST distribution is possible. Next NOPN value:

end do outer

call mma_deallocate(ISM)

contains

subroutine ErrorTrap()

  integer(kind=iwp) NOPN, ISYM

  write(u6,*) ' MKCNF ERROR: Unforeseen calamity.'
  write(u6,*) ' At end of loop over NOPN, the number of'
  write(u6,*) ' configurations generated does not match'
  write(u6,*) ' that which was allocated.'
  write(u6,*) ' INFO table in ICNFTAB says:'
  write(u6,*)
  write(u6,*) '  NOPN ISYM       Nr of conf Start point  Words/config'
  do NOPN=MINOP,MAXOP
    NCLS = (NEL-NOPN)/2
    NOCC = NCLS+NOPN
    do ISYM=1,NSYM
      NCNF = ICNFTAB(KINFO+0+3*(ISYM-1+NSYM*(NOPN-MINOP)))
      KCNFSTA = ICNFTAB(KINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))
      LENCNF = ICNFTAB(KINFO+2+3*(ISYM-1+NSYM*(NOPN-MINOP)))
      write(u6,'(1X,2I4,5X,3I12)') NOPN,ISYM,NCNF,KCNFSTA,LENCNF
    end do
  end do
  call ABEND()

end subroutine ErrorTrap

end subroutine MKCONF
