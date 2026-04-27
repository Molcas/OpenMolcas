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

subroutine SYGTOSD(ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,CISYG,CISD,detocc,detcoeff,SPTRA)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICNFTAB(*), ISPNTAB(*), ISSTAB(*), IFSBTAB(*)
real(kind=wp), intent(in) :: CISYG(*), SPTRA(*)
real(kind=wp), intent(_OUT_) :: CISD(*), detcoeff(*)
character(len=*), intent(_OUT_) :: detocc(*)
integer(kind=iwp) :: I, IBLK, ICNF, idet, IEL, IEL1, IEL2, IERR, IFORM, IFSB, IMORS, IOCC, IOEND, IORB, IOSTA, IPART, IPOS, IREST, &
                     ISBSTR, ISORB, ISPART, ISPD, ISPEND, ISPN, ISPSTA, ISST, ISUM, ISYGEND, ISYGSTA, IWORD, IWRD, JSST, KCNF, &
                     KCNFINF, KFSB, KGSLIM, KGSORB, KHSHMAP, KMRSSBS, KSPN, KSPNINF, KSSTARR, KSSTTB, LSPTRA, LSYM, MAXOP, MINOP, &
                     MORSBITS, MXBLK, NACTEL, NAPART, NASPRT, NBLK, NCLSD, NCNF, NCPL, NFSB, NHEAD, NHSHMAP, NO, NOCC, NOP, NOPEN, &
                     NORB, NSP, NSPD, NSSTP, NSYM, NWRD
integer(kind=iwp), allocatable :: ndim(:), OccArr(:), OrbArr(:), SBSET(:), SSArr(:), STArr(:)
real(kind=wp), allocatable :: BLK(:)
character, allocatable :: occ(:)
integer(kind=iwp), external :: OCC2MRS

! Unbutton the configuration table:
NACTEL = ICNFTAB(3)
NORB = ICNFTAB(4)
MINOP = ICNFTAB(5)
MAXOP = ICNFTAB(6)
NSYM = ICNFTAB(7)
LSYM = ICNFTAB(8)
NAPART = ICNFTAB(9)
IFORM = ICNFTAB(10)
NHEAD = 10
KGSORB = NHEAD+1
KGSLIM = KGSORB+(NSYM+1)*(NAPART+1)
KCNFINF = KGSLIM+2*NAPART
!TEST write(u6,*) ' Table at KGSORB:'
!TEST do ipart=0,napart
!TEST   write(u6,'(1x,10i5)') (icnftab(kgsorb+isym+(nsym+1)*ipart),isym=0,nsym)
!TEST end do
!TEST write(u6,*) ' Table at KGSLIM:'
!TEST write(u6,'(1x,10i5)') (icnftab(kgslim+2*(ipart-1)),ipart=1,napart)
!TEST write(u6,'(1x,10i5)') (icnftab(kgslim+1+2*(ipart-1)),ipart=1,napart)
! Unbutton the spin coupling table:
KSPNINF = 9
! Unbutton the Substring Table:
NASPRT = ISSTAB(5)
MORSBITS = ISSTAB(6)
NSSTP = ISSTAB(7)
KSSTTB = 15
!TEST KSBSMRS = ISSTAB(11)
KMRSSBS = ISSTAB(12)
! Unbutton the Fock Sector Block table:
NHEAD = 7
KSSTARR = NHEAD+1
NFSB = IFSBTAB(2)
NHSHMAP = IFSBTAB(6)
KHSHMAP = IFSBTAB(7)

! MXBLK=Largest individual SYG block of determinants:
MXBLK = 0
idet = 0
do NOPEN=MINOP,MAXOP
  NCNF = ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP)))
  NCPL = ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+2)
  MXBLK = max(NCNF*NCPL,MXBLK)
end do
! A variety of small temporary arrays used:
call mma_allocate(BLK,MXBLK,Label='BLK')
call mma_allocate(ORBARR,NACTEL,Label='OrbArr')
call mma_allocate(OCCARR,2*NORB,Label='OccArr')
call mma_allocate(STARR,NASPRT,Label='STArr')
call mma_allocate(NDIM,NASPRT,Label='nDim')
call mma_allocate(SSARR,NASPRT,Label='SSArr')
call mma_allocate(SBSET,NSSTP,Label='SBSet')
call mma_allocate(occ,norb,label='occ')
! We will need later the accumulated number of substrings of
! earlier substring types:
ISUM = 0
do ISST=1,NSSTP
  SBSET(ISST) = ISUM
  ISUM = ISUM+ISSTAB(KSSTTB+5*(ISST-1))
end do
! Loop over nr of open shells.
ISYGEND = 0
do NOPEN=MINOP,MAXOP
  NCLSD = (NACTEL-NOPEN)/2
  if (NCLSD < 0) cycle
  if (2*NCLSD+NOPEN /= NACTEL) cycle
  NOCC = NCLSD+NOPEN
  if (NOCC > NORB) cycle
  NCNF = ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP)))
  if (NCNF == 0) cycle
  KCNF = ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+1)
  NWRD = ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+2)
  NCPL = ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+1)
  NSPD = ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+2)
  KSPN = ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+4)
  if (NSPD == 0) cycle
  if (NCPL == 0) cycle
  ! ISYGSTA=1st element of each block
  NBLK = NCPL*NCNF
  ISYGSTA = ISYGEND+1
  ISYGEND = ISYGEND+NBLK
  ! Location of spin coupling coefficients:
  LSPTRA = ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+5)
  ! Matrix multiplication into temporary array:
  call DGEMM_('N','N',NSPD,NCNF,NCPL,One,SPTRA(LSPTRA),NSPD,CISYG(ISYGSTA),NCPL,Zero,BLK,NSPD)

  ! There is no phase factor in the reorder of orbitals from SYG
  ! to SD. But there is a fairly lengthy procedure for finding the
  ! correct position of each CI coefficient.
  IBLK = 0
  ! Loop over configurations
  IWORD = 0 ! dummy initialize
  do ICNF=1,NCNF
    if (IFORM == 1) then
      ORBARR(1:NOCC) = ICNFTAB(KCNF+NWRD*(ICNF-1):KCNF+NWRD*(ICNF-1)+NOCC-1)
    else if (IFORM == 2) then
      IEL2 = 0
      IEL1 = NCLSD
      do IORB=1,NORB
        IOCC = ICNFTAB(KCNF-1+IORB+NWRD*(ICNF-1))
        if (IOCC == 1) then
          IEL1 = IEL1+1
          ORBARR(IEL1) = IORB
        else
          IEL2 = IEL2+1
          ORBARR(IEL2) = IORB
        end if
      end do
    else if (IFORM == 3) then
      do IEL=1,NOCC
        IWRD = (3+IEL)/4
        IREST = (3+IEL)-4*IWRD
        if (IREST == 0) IWORD = ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
        IORB = mod(IWORD,256)
        IWORD = IWORD/256
        ORBARR(IEL) = IORB
      end do
    else if (IFORM == 4) then
      IEL2 = 0
      IEL1 = NCLSD
      do IORB=1,NORB
        IWRD = (IORB+14)/15
        IREST = IORB+14-15*IWRD
        if (IREST == 0) IWORD = ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
        IOCC = mod(IWORD,4)
        IWORD = IWORD/4
        if (IOCC == 1) then
          IEL1 = IEL1+1
          ORBARR(IEL1) = IORB
        else
          IEL2 = IEL2+1
          ORBARR(IEL2) = IORB
        end if
      end do
    end if

    !TEST write(u6,'(1x,a,10i5)') 'Configuration:',(ORBARR(iel),iel=1,nocc)
    !TEST write(u6,*) ' Loop over spin determinants.'
    ! Loop over spin determinants
    do ISPD=1,NSPD
      !TEST write(u6,'(1x,a,20i3)') 'Spin determinant:',(ISPNTAB(KSPN-1+I+NOPEN*(ISPD-1)),I=1,nopen)
      IBLK = IBLK+1
      ! count the determinants
      idet = idet+1
      ! Construct occupation number array:
      OCCARR(:) = 0
      do IEL=1,NCLSD
        IORB = ORBARR(IEL)
        OCCARR(2*IORB-1) = 1
        OCCARR(2*IORB) = 1
      end do
      do I=1,NOPEN
        ! Spin of each electron is coded as 1 for alpha, 0 for beta.
        ISPN = ISPNTAB(KSPN-1+I+NOPEN*(ISPD-1))
        IEL = NCLSD+I
        ISORB = 2*ORBARR(IEL)-ISPN
        OCCARR(ISORB) = 1
      end do
      ! Identify substrings:
      ! Loop over active partitions. Subdivide as needed into subpartitions.
      !TEST write(u6,*) ' Identify substrings.'
      !TEST write(u6,'(1x,a,10i5)') 'Occupation array:',(OCCARR(isorb),isorb=1,2*norb)
      ! construct occupation array in 0,u,d,2 format
      do IORB=1,2*norb-1,2
        if ((OCCARR(IORB) == 1) .and. (OCCARR(1+IORB) == 1)) then
          occ((IORB+1)/2) = '2'
        else if ((OCCARR(IORB) == 1) .and. (OCCARR(1+IORB) == 0)) then
          occ((IORB+1)/2) = 'u'
        else if ((OCCARR(IORB) == 0) .and. (OCCARR(1+IORB) == 1)) then
          occ((IORB+1)/2) = 'd'
        else if ((OCCARR(IORB) == 0) .and. (OCCARR(1+IORB) == 0)) then
          occ((IORB+1)/2) = '0'
        end if
      end do

      IOEND = 0
      ISPEND = 0
      !TEST write(u6,'(1x,a,10i5)') 'NAPART:',NAPART
      do IPART=1,NAPART
        !TEST write(u6,'(1x,a,10i5)') ' In loop, IPART:',IPART
        NOP = 2*ICNFTAB(KGSORB+(NSYM+1)*IPART)
        !TEST write(u6,'(1x,a,10i5)') '            NOP:',NOP
        !TEST if (NOP == 0) write(u6,*) ' (Skip it.)'
        if (NOP == 0) cycle
        NSP = (NOP+MORSBITS-1)/MORSBITS
        !TEST write(u6,'(1x,a,10i5)') '            NSP:',NSP
        ISPSTA = ISPEND+1
        ISPEND = ISPEND+NSP
        do ISPART=ISPSTA,ISPEND
          !TEST write(u6,'(1x,a,10i5)') 'Loop lims ISPSTA,ISPEND:',ISPSTA,ISPEND
          NO = min(NOP,MORSBITS)
          NOP = NOP-NO
          IOSTA = IOEND+1
          IOEND = IOEND+NO
          !TEST write(u6,'(1x,a,10i5)') 'IOSTA,IOEND:',IOSTA,IOEND
          !TEST write(u6,'(1x,a,10i5)') 'Occ array:',(OCCARR(ISORB),ISORB=IOSTA,IOEND)
          IMORS = OCC2MRS(NO,OCCARR(IOSTA))
          !TEST write(u6,'(1x,a,10i5)') 'IMORS=',IMORS
          ! Position in Morsel-to-Substring table:
          IPOS = KMRSSBS+2*(IMORS+(2**MORSBITS)*(ISPART-1))
          ! Substring ID number
          ISBSTR = ISSTAB(IPOS)
          ! Test:
          !TEST JMORS = ISSTAB(KSBSMRS+2*(ISBSTR-1))
          !TEST if (IMORS /= JMORS) then
          !TEST   write(u6,*) ' Mistranslated morsel!!'
          !TEST   write(u6,'(1x,a,4i12)') 'IMORS->ISBSTR:',IMORS,ISBSTR
          !TEST   write(u6,'(1x,a,4i12)') 'but ISBSTR->IMORS:',ISBSTR,JMORS
          !TEST   write(u6,'(1x,a,4i12)') 'KMRSSBS:',KMRSSBS
          !TEST   write(u6,'(1x,a,4i12)') 'KSBSMRS:',KSBSMRS
          !TEST   call ABEND()
          !TEST end if
          !TEST write(u6,'(1x,a,10i5)') 'ISBSTR:',ISBSTR
          ! Substring type ISST, nr of such substrings is NDIM
          ISST = ISSTAB(IPOS+1)
          !TEST write(u6,'(1x,a,10i5)') 'ISST  :',ISST
          STARR(ISPART) = ISST
          ndim(ISPART) = ISSTAB(KSSTTB+5*(ISST-1))
          SSARR(ISPART) = ISBSTR-SBSET(ISST)
        end do
      end do
      !TEST write(u6,*) ' Finally, substring types and substrings:'
      !TEST write(u6,'(1x,a,10i5)') 'Substr types:',(STARR(ISPART),ISPART=1,NASPRT)
      !TEST write(u6,'(1x,a,10i5)') 'Substrings  :',(SSARR(ISPART),ISPART=1,NASPRT)
      !TEST write(u6,'(1x,a,10i5)') 'Dimensions  :',(NDIM(ISPART),ISPART=1,NASPRT)
      ! Position within FS block:
      IPOS = (SSARR(NASPRT)-1)
      do ISPART=NASPRT-1,1,-1
        IPOS = ndim(ISPART)*IPOS+(SSARR(ISPART)-1)
      end do
      IPOS = IPOS+1
      ! Identify Fock Sector Block:
      !TEST write(u6,*) ' Arguments in HSHGET call:'
      !TEST write(u6,'(1x,a,10i5)') 'Key:',(STARR(ISPART),ISPART=1,NASPRT)
      !TEST write(u6,'(1x,a,10i5)') 'Size of key:',NASPRT
      !TEST write(u6,'(1x,a,10i5)') 'Size of items stored:',NASPRT+2
      !TEST write(u6,'(1x,a,10i5)') 'Items stored at KSSTARR=',KSSTARR
      !TEST write(u6,'(1x,a,10i5)') '      Map size  NHSHMAP=',NHSHMAP
      !TEST write(u6,'(1x,a,10i5)') '  Map stored at KHSHMAP=',KHSHMAP
      call HSHGET(STARR,NASPRT,NASPRT+2,IFSBTAB(KSSTARR),NHSHMAP,IFSBTAB(KHSHMAP),IFSB)
      !TEST write(u6,'(1x,a,10i5)') ' Map returns index IFSB=',IFSB
      !TEST write(u6,'(1x,a,10i5)') ' Item stored there is  =',(IFSBTAB(KSSTARR-1+ISPART+(NASPRT+2)*(IFSB-1)),ISPART=1,NASPRT+2)
      ! Position of this FS block in SD wave function:
      KFSB = IFSBTAB(KSSTARR+(NASPRT+2)*IFSB-1)
      ! Temporary check, may be removed later. See that we have picked up
      ! the correct FS block.
      IERR = 0
      do ISPART=1,NASPRT
        JSST = IFSBTAB(KSSTARR-1+ISPART+(NASPRT+2)*(IFSB-1))
        ISST = STARR(ISPART)
        if (ISST /= JSST) IERR = 1
      end do
      if (IERR /= 0) then
        write(u6,*) ' SYGTOSD Error: Hash map returned the wrong FS block!'
        write(u6,'(1x,a,8I8)') 'NOPEN,ICNF,ISPN:',NOPEN,ICNF,ISPN
        write(u6,'(1x,a,20I3)') 'Configuration:',(ORBARR(IEL),IEL=1,NACTEL)
        write(u6,'(1x,a,20I3)') 'Determinant:',(OCCARR(ISORB),ISORB=1,2*NORB)
        write(u6,'(1x,a,10I5)') 'Substring type combination:',(STARR(ISPART),ISPART=1,NASPRT)
        write(u6,'(1x,a,10I5)') 'Substring combination:',(SSARR(ISPART),ISPART=1,NASPRT)
        write(u6,'(1x,a,8I8)') 'Hash table says IFSB=',IFSB
        if ((IFSB > 0) .and. (IFSB <= NFSB)) then
          write(u6,'(1x,a,8I8)') 'but that FS block would contain',(IFSBTAB(KSSTARR-1+ISPART+(NASPRT+2)*(IFSB-1)),ISPART=1,NASPRT)
        else
          write(u6,*) 'but there is no such FS block!'
          write(u6,*) ' The FS block table follows:'
          call PRFSBTAB(IFSBTAB)
        end if
        call ABEND()
      end if
      ! Finally:
      CISD(KFSB-1+IPOS) = BLK(IBLK)
      detcoeff(idet) = BLK(IBLK)
      write(detocc(idet),*) occ
      ! End of spin-determinant loop
    end do
    ! End of loop over configurations
  end do
  ! End of loop over nr of open shells
end do
call mma_deallocate(BLK)
call mma_deallocate(ORBARR)
call mma_deallocate(OCCARR)
call mma_deallocate(STARR)
call mma_deallocate(NDIM)
call mma_deallocate(SSARR)
call mma_deallocate(SBSET)
call mma_deallocate(occ)

return

end subroutine SYGTOSD
