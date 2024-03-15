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

subroutine get_Umn(PHP,EnIn,DHAM,IPCSF,IPCNF,MXPDIM,DTOC,IPRODT,ICONF,IREFSM,ONEBOD,ECORE,NACTOB,NCONF,NEL,NAEL,NBEL,NPCSF,NPCNF, &
                   TUVX,iterSplit,ITER,NTEST,ExFac,IREOTS)
! ARGUMENTS :
! ===========
! PHP    : AA Block Hamiltonian un-dressed                  (Output)
! EnIn   : energy value of the root selected                (Input)
! DHAM   : Dressed AA block Hamiltonian                     (Output)
! IPCSF  : CSF's order - Index Array -                      (Input)
! IPCNF  : CNF's order - Index Array -                      (Input)
! MXPDIM : Total number of CSFs                             (Input)
! DTOC   : Transformation matrix between CSF's and DET's    (Input)
! IPRODT : Prototype determinants                           (Input)
! ICONF  : List of configurations                           (Input)
! IREFSM : symmetry of considered CI space                  (Input)
! Onebod : one body hamilton matrix in rectangular form     (Input)
! ECORE  : Core energy                                      (Input)
! NACTOB : Number of active orbitals                        (Input)
! NCONF  : Number of CNFs of symmetry IREFSM                (Input)
! NEL    : total number of active electrons                 (Input)
! NAEL   : number of alpha active electron                  (Input)
! NBEL   : number of beta active electron                   (Input)
! NPCSF  : Number of CSFs in AA block                       (Input)
! NPCNF  : Number of CNFs in AA block                       (Input)
! TUVX   : Two-electron integrals (MO space)                (Input)
! NTEST  :
! ExFac  :
! IREOTS : Type => symmetry reordering array

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MXPDIM, NCONF, IPCSF(MXPDIM), IPCNF(NCONF), IPRODT(*), ICONF(*), IREFSM, NACTOB, NEL, NAEL, NBEL, &
                                 NPCSF, NPCNF, iterSplit, ITER, IREOTS(NACTOB)
real(kind=wp), intent(in) :: EnIn, DTOC(*), ONEBOD(NACTOB,NACTOB), ECORE, TUVX(*), ExFac
real(kind=wp), intent(out) :: PHP(NPCSF*(NPCSF+1)/2), DHAM(NPCSF*(NPCSF+1)/2)
integer(kind=iwp), intent(inout) :: NTEST
integer(kind=iwp) :: iAlpha, IATYP, IIA, IIAB, IIL, IILACT, IILB, IIR, IIRACT, IIRB, IIRMAX, ILAI, &
                     ILRI, ILRO, ILTYP, IRTYP, ITYP, Mindex, MXCSFC, MXXWS, NCSFA, &
                     NCSFL, NCSFR, Nindex
real(kind=wp), allocatable :: AuxC(:,:), AuxD(:), AuxV(:,:), Scr(:)
integer(kind=iwp), allocatable:: ICNL(:), ICNR(:), ICNQ(:)
real(kind=wp), allocatable:: CNHCNM(:), PHPS(:)
#include "spinfo.fh"

if (NTEST >= 30) then
  write(u6,*) ' Input in get_Umn'
  write(u6,*) ' =================='
  write(u6,*) ' Number of CNFs ',NCONF
  write(u6,*) ' Number of CSFs ',MXPDIM
  write(u6,*) ' Configurations included :'
  call IWRTMA(IPCNF,1,NCONF,1,NCONF)
  write(u6,*) ' CSFs included :'
  call IWRTMA(IPCSF,1,MXPDIM,1,MXPDIM)
  write(u6,*) ' Number of CNFs in AA block:',NPCNF
  write(u6,*) ' Number of CSFs in AA block:',NPCSF
end if

DHAM(:) = Zero
! construct the Dressed Hamiltonian matrix
MXCSFC = 0
do ITYP=1,NTYP
  MXCSFC = max(MXCSFC,NCSFTP(ITYP))
end do
!write(u6,*) 'MXCSFC = ',MXCSFC
call mma_allocate(AuxD,MXCSFC,label='AuxDia')
call mma_allocate(AuxV,NPCSF,MXCSFC,label='AuxVer')
call mma_allocate(AuxC,NPCSF,MXCSFC,label='AuxCopy')

Call mma_allocate(ICNL,NEL,Label='ICNL')
Call mma_allocate(ICNR,NEL,Label='ICNR')
Call mma_allocate(ICNQ,NEL,Label='ICNQ')

Call mma_allocate(CNHCNM,MXCSFC**2,Label='CNHCNM')
Call mma_allocate(PHPS,MXCSFC**2,Label='PHPS')

call mma_maxDBLE(MXXWS)
MXXWS = MXXWS/2
call mma_allocate(Scr,MXXWS,label='EXHSCR')


! Shape of PHP matrix:
!
!           CNFL (n)
!
!      1  2  4  7 11 16 22 ...
!  C   -  3  5  8 12 17 23 ...
!  N   -  -  6  9 13 18 24 ...
!  F   -  -  - 10 14 19 25 ...
!  R   -  -  -  - 15 20 26 ...
! (m)  -  -  -  -  - 21 27 ...
!      -  -  -  -  -  - 28 ...
!      -  -  -  -  -  -  - ...

!***********************************************************************
! Now a loop over alpha will start to calculate:                       *
!   1) BB-Block DIAGONAL element 1/(En-H(alpha,alpha))                 *
!   2) AB-Block array H(m,alpha)
!***********************************************************************
if ((ITER /= 1) .or. (iterSplit /= 1)) then
  IIAB = 1
  do iAlpha=NPCNF+1,NCONF ! Loop over alpha
    CNHCNM(:)=0.0D0
    !write(u6,*) 'iAlpha = ',iAlpha
    call GETCNF_LUCIA(ICNL,IATYP,IPCNF(iAlpha),ICONF,IREFSM,NEL)
    NCSFA = NCSFTP(IATYP)
    !write(u6,*) 'NCSFA = ',NCSFA
    call CNHCN(ICNL,IATYP,ICNL,IATYP,CNHCNM,SCR,NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX, &
               NTEST,ExFac,IREOTS)
    do IIA=1,NCSFA
      ILAI = IIA*IIA
      if (NTEST >= 30) write(u6,*) 'ILAI =',ILAI
      !AuxD(IIA) = CNHCNM(ILAI)
      !write(u6,*) 'AuxD(IIA)',AuxD(IIA)
      AuxD(IIA) = One/(EnIn-CNHCNM(ILAI))
      if (NTEST >= 30) write(u6,*) 'AuxD(IIA)',AuxD(IIA)
    end do
    !*************** 2) AB-Block Array (alpha Column) ********************
    IILB = 1
    do Mindex=1,NPCNF ! Loop over AB-Block
      CNHCNM(:)=0.0D0
      !write(u6,*) 'Mindex',Mindex
      call GETCNF_LUCIA(ICNR,ILTYP,IPCNF(Mindex),ICONF,IREFSM,NEL)
      NCSFL = NCSFTP(ILTYP)
      !write(u6,*) 'NCSFL = ',NCSFL
      call CNHCN(ICNL,IATYP,ICNR,ILTYP,CNHCNM,SCR,NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB, &
                 TUVX,NTEST,ExFac,IREOTS)
      if (NTEST >= 30) then
        write(u6,*) 'M_Alpha elements'
        call wrtmat(CNHCNM,MXCSFC,MXCSFC,MXCSFC,MXCSFC)
      end if
      do IIL=1,NCSFL
        do IIA=1,NCSFA
          IILACT = IILB-1+IIL
          ILAI = (IIL-1)*NCSFA+IIA
          !ILAI = (IIA-1)*MXCSFC+IIL
          !ILAI = (IIL-1)*MXCSFC+IIA
          AuxV(IILACT,IIA) = CNHCNM(ILAI)
          !write(u6,*) 'ILAI, ILACT, IIA = ',ILAI,ILACT,IIA
          if (NTEST >= 30) write(u6,*) 'AuxV(IILACT,IIA)',AuxV(IILACT,IIA)
          AuxC(IILACT,IIA) = AuxV(IILACT,IIA)*AuxD(IIA)
          if (NTEST >= 30) write(u6,*) 'AuxC(IILACT,IIA)',AuxC(IILACT,IIA)
        end do
      end do
      IILB = IILB+NCSFL
    end do ! End loop over AB-Block
    if (NTEST >= 30) then
      write(u6,*) 'AB-Block Vertical Vector'
      call wrtmat(AuxV,NPCSF,NCSFA,NPCSF,NCSFA)
      write(u6,*) 'AB-Block Vertical Vector times Daa'
      call wrtmat(AuxC,NPCSF,NCSFA,NPCSF,NCSFA)
    end if
    !*********************************************************************
    call dGeMM_Tri('N','T',NPCSF,NPCSF,NCSFA,One,AuxC,NPCSF,AuxV,NPCSF,One,DHAM,NPCSF)
    if (NTEST >= 30) call TRIPRT('correction to the AA block',' ',DHAM,NPCSF)
    IIAB = IIAB+NCSFA
  end do ! End of the loop over iAlpha
end if
!***********************************************************************
!    1. A-Block matrix element                                         *
!***********************************************************************
IILB = 1
do Nindex=1,NPCNF ! Loop over the AA-block (vertical index)
  if (NTEST >= 30) write(u6,*) 'Nindex',Nindex
  !write(u6,*) 'IILB',IILB
  call GETCNF_LUCIA(ICNR,ILTYP,IPCNF(Nindex),ICONF,IREFSM,NEL)
  NCSFL = NCSFTP(ILTYP)
  !write(u6,*) 'NCSFL = ',NCSFL

  IIRB = 1
  do Mindex=1,Nindex ! Loop over the AA-block (horizontal index)
    call FZero(PHPS,MXCSFC*MXCSFC)
    !write(u6,*) 'Nindex,Mindex',Nindex,Mindex
    call GETCNF_LUCIA(ICNQ,IRTYP,IPCNF(Mindex),ICONF,IREFSM,NEL)
    NCSFR = NCSFTP(IRTYP)
    !write(u6,*) 'NCSFR = ',NCSFR
    call CNHCN(ICNR,ILTYP,ICNQ,IRTYP,PHPS,SCR,NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX, &
               NTEST,ExFac,IREOTS)
    if (NTEST >= 30) then
      write(u6,*) 'AA block elements'
      call wrtmat(PHPS,MXCSFC,MXCSFC,MXCSFC,MXCSFC)
    end if
    do IIL=1,NCSFL
      if (IILB == IIRB) then
        IIRMAX = IIL
      else
        IIRMAX = NCSFR
      end if
      do IIR=1,IIRMAX
        IIRACT = IIRB-1+IIR
        IILACT = IILB-1+IIL
        !ILRI = (IIR-1)*MXCSFC+IIL
        ILRI = (IIR-1)*NCSFL+IIL
        !^ Forse questo e' quello giusto; la precedente formula
        !  potrebbe essere fonte di BUGS! Vedremo!
        ILRO = ((IILACT*IILACT-IILACT)/2)+IIRACT
        PHP(ILRO) = PHPS(ILRI)
        !write(u6,*) 'ILRI, ILRO = ',ILRI,ILRO
        !write(u6,*) 'PHP(ILRO)',PHP(ILRO)
      end do
    end do
    IIRB = IIRB+NCSFR
  end do ! End loop over the AA-block (horizontal index)
  IILB = IILB+NCSFL
end do ! End loop over the AA-block (vertical index)

!***********************************************************************
! Let's add Hmn (PHP) to the correction (DHAM)
DHAM(:) = DHAM(:)+PHP(:)
!***********************************************************************
if (NTEST >= 30) then
  write(u6,*) 'AA-Block matrix un-dressed'
  call wrtmat(PHP,NPCSF*(NPCSF+1)/2,1,NPCSF*(NPCSF+1)/2,1)
  write(u6,*) 'AA-Block matrix dressed'
  call wrtmat(DHAM,NPCSF*(NPCSF+1)/2,1,NPCSF*(NPCSF+1)/2,1)
end if

if (NTEST >= 30) then
  call TRIPRT('AA block Hamiltonian Matrix un-dressed',' ',PHP,NPCSF)
  call TRIPRT('Dressed AA block Hamiltonian Matrix',' ',DHAM,NPCSF)
end if
call mma_deallocate(AuxD)
call mma_deallocate(AuxV)
call mma_deallocate(AuxC)
call mma_deallocate(ICNL)
call mma_deallocate(ICNR)
call mma_deallocate(ICNQ)
call mma_deallocate(CNHCNM)
call mma_deallocate(PHPS)
call mma_deallocate(Scr)

end subroutine get_Umn
