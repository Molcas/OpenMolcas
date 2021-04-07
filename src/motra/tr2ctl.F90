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

subroutine TR2CTL(CMO)
! Two-electron integral transformation program: control section
!
! Purpose: Set up of memory locations (decide if out of core is
!          needed. Loop over symmetry blocks of AO integrals.
!          The transformation routine TRAMO is called for each
!          symmetry block of integrals.

implicit real*8(A-H,O-Z)
#include "motra_global.fh"
#include "trafo_motra.fh"
#include "files_motra.fh"
#include "WrkSpc.fh"
real*8 CMO(*)
dimension NBSX(8), KEEP(8), ISTSQ(8)
logical FoundTwoEls, DoDirect, DoCholesky, ISQUAR

! Set time at start of transformation

call SETTIM()

! Initiate unit LUTWOMO (two-electron integrals in MO basis)

call DANAME_MF(LUTWOMO,FNTWOMO)
IAD13 = 0
call iCopy(nTraToc,[0],0,iTraToc,1)
call iDAFILE(LUTWOMO,1,iTraToc,nTraToc,IAD13)

! Initiate unit LUTWOAO (two-electron integrals in AO basis)

call f_Inquire(FnTwoAo,FoundTwoEls)
call DecideOnDirect(.false.,FoundTwoEls,DoDirect,DoCholesky)
if (.not. DoCholesky) call OPNORD(IRC,0,FNTWOAO,LUTWOAO)

call GETORD(IRC,ISQUAR,NSYM2,NBSX,KEEP)

! Compare content of 1el and 2el integral file

if (NSYM2 /= NSYM) then
  write(6,*) 'Tr2Ctl: NSYM2.NE.NSYM'
  write(6,*) 'NSYM2=',NSYM2
  write(6,*) 'NSYM=',NSYM
  call Abend()
end if
do ISYM=1,NSYM
  NB1 = NBAS(ISYM)
  NB2 = NBSX(ISYM)
  if (NB1 /= NB2) then
    write(6,*) 'Tr2Ctl: NB1.NE.NB2'
    write(6,*) 'NB1=',NB1
    write(6,*) 'NB2=',NB2
    call Abend()
  end if
end do

! Precompute start points

ISTBS = 1
do ISYM=1,NSYM
  ISTSQ(ISYM) = ISTBS+NBAS(ISYM)*NFRO(ISYM)
  ISTBS = ISTBS+NBAS(ISYM)*NBAS(ISYM)
end do

! Loop over quadruples of symmetries (nsp,nsq,nsr,nss)
! Note that the integrals on LUTWOAO have to be sorted in the
! same order as the loop structure below.

if (iPrint >= 0) write(6,2000)
2000 format(/7X,'SYMMETRY',2X,'BASIS FUNCTIONS',6X,' ORBITALS',6X,'INTEGRALS   CPU(SEC)  I/O(SEC)')
IBATCH = 0
do NSP=1,NSYM
  NBP = NBAS(NSP)
  NOP = NORB(NSP)
  KEEPP = KEEP(NSP)
  LMOP = ISTSQ(NSP)
  ISP = NSP
  do NSQ=1,NSP
    NBQ = NBAS(NSQ)
    NOQ = NORB(NSQ)
    KEEPQ = KEEP(NSQ)
    LMOQ = ISTSQ(NSQ)
    NSPQ = ieor(NSP-1,NSQ-1)+1
    ISQ = NSQ
    do NSR=1,NSP
      NBR = NBAS(NSR)
      NOR = NORB(NSR)
      KEEPR = KEEP(NSR)
      LMOR = ISTSQ(NSR)
      NSPQR = ieor(NSPQ-1,NSR-1)+1
      ISR = NSR
      NSSM = NSR
      if (NSR == NSP) NSSM = NSQ
      do NSS=1,NSSM
        NBS = NBAS(NSS)
        NOS = NORB(NSS)
        KEEPS = KEEP(NSS)
        LMOS = ISTSQ(NSS)
        NSPQRS = ieor(NSPQR-1,NSS-1)+1
        if (NSPQRS /= 1) goto 101
        IBATCH = IBATCH+1
        ISS = NSS

        ! Check the loop conditions and skip transformation step if possible

        KEEPT = KEEPP+KEEPQ+KEEPR+KEEPS
        NORBP = NOP*NOQ*NOR*NOS
        if (NORBP == 0) goto 101
        if (KEEPT /= 0) then
          write(6,*) 'Tr2Ctl: NORBP /= 0 .AND. KEEPT /= 0'
          write(6,*) 'NORBP=',NORBP
          write(6,*) 'KEEPT=',KEEPT
          call Abend()
        end if

        ! Define matrix sizes

        if (ISR == ISS) then
          NBPQ = (NBP+NBP**2)/2
          NBRS = (NBR+NBR**2)/2
          NOVX = (NOR+NOR**2)/2
        else
          NBPQ = NBP*NBQ
          NBRS = NBR*NBS
          NOVX = NOR*NOS
        end if

        ! Set up dynamic memory

        ! INTBUF = Size of output buffer.
        !INTBUF = 256*256
        INTBUF = 16*256*256
        NW1 = 2*KBUF
        NW2 = max(INTBUF,NBP*NOQ,NBQ*NOP)
        ! NW2 is size of 'X1' in TRAMO, used as LBUF in call to RDORD and RDORD_.
        ! LBUF-1 must be at least 'klB':
        NW2 = max(NW2,NBRS+1,NBPQ+1)
        NW3 = max(NBR**2,NBP**2,NOQ*NBP,NOVX)
        NW4 = max(NBR*NOS,NBQ*NOP)
        call GETMEM('OUTBUF','ALLO','REAL',LW1,NW1)
        call GETMEM('X1','ALLO','REAL',LW2,NW2)
        call GETMEM('X2','ALLO','REAL',LW3,NW3)
        call GETMEM('X3','ALLO','REAL',LW4,NW4)
        call GetMem('iDsk','Allo','Inte',ipiDsk,3*nOVX)
        call GETMEM('VXPQ','MAX','REAL',LW5,MEMX)

        if (DoCholesky) then ! save some space for GenInt
          MEMX = max(MEMX-MEMX/10,0)
        end if

        call GETMEM('VXPQ','ALLO','REAL',LW5,MEMX)

        if (MEMX < NOVX) then
          write(6,*) 'Tr2Ctl: MEMX.LT.NOVX'
          write(6,*) 'MEMX=',MEMX
          write(6,*) 'NOVX=',NOVX
          call Abend()
        end if

        ! transform the symmetry block (ISP,ISQ|ISR,ISS)

        iTraToc(IBATCH) = IAD13
        ! NW2 is size of 'X1' in TRAMO, used as LBUF in call to RDORD and RDORD_.
        call TRAMO(NW2,WORK(LW1),nW1,WORK(LW2),nW2,WORK(LW3),nW3,WORK(LW4),nW4,WORK(LW5),MEMX,CMO,iWork(ipiDsk),nOVX)
        call TIMING(CPT,CPE,TIOT,TIOE)
        if (iPrint >= 0) write(6,2100) ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,NOP,NOQ,NOR,NOS,LTUVX,CPE,TIOE
2100    format(7X,4I2,1X,4I4,2X,4I4,3X,I9,F11.2,F10.2)
        call Xflush(6)

        ! Deallocate work space

        call GETMEM('VXPQ','FREE','REAL',LW5,MEMX)
        call GetMem('iDsk','Free','Inte',ipiDsk,3*nOVX)
        call GETMEM('X3','FREE','REAL',LW4,NW4)
        call GETMEM('X2','FREE','REAL',LW3,NW3)
        call GETMEM('X1','FREE','REAL',LW2,NW2)
        call GETMEM('OUTBUF','FREE','REAL',LW1,NW1)

        ! End of loop over quadruples of symmetries

101     continue
      end do
    end do
  end do
end do
call TIMING(CPT,CPE,TIOT,TIOE)
if (iPrint >= 0) write(6,2200) CPT,TIOT
2200 format(/6X,' TOTAL CPU TIME(SEC)',F8.2,'TOTAL I/O TIME(SEC)',F8.2)

! Close LUTWOAO

if (.not. DoCholesky) then
  call CLSORD(IRC,0)
end if

! Write address block to LUTWOMO

IAD13 = 0
call iDAFILE(LUTWOMO,1,iTraToc,nTraToc,IAD13)
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(6,'(6X,A)') 'DISK ADRESSES FOR SYMMETRY BLOCKS OF TRANSFORMED TWO-ELECTRON INTEGRALS'
  write(6,'(6X,10I8)') (iTraToc(I),I=1,nTraToc)
end if

call DACLOS(LUTWOMO)

return

end subroutine TR2CTL
