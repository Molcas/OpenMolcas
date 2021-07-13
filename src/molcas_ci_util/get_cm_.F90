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

subroutine get_Cm_(IPCSF,IPCNF,MXPDIM,NCONF,NPCSF,NPCNF,Cn,lrootSplit,EnFin,DTOC,IPRODT,ICONF,IREFSM,ONEBOD,ECORE,NACTOB,NEL,NAEL, &
                   NBEL,DIAG,TUVX,NTEST,ExFac,IREOTS,Ctot)
! Obtain Cm coefficients out of the AA Block for root = lrootSplit
! to the first order in Lowdin equation
!
! ARGUMENTS :
! ===========
! IPCSF      : CSF's order - Index Array -                    (Input)
! IPCNF      : CNF's order - Index Array -                    (Input)
! MXPDIM     : Total number of CSFs                           (Input)
! NCONF      : Total Number of CNFs of symmetry IREFSM        (Input)
! NPCSF      : Number of CSFs in AA block                     (Input)
! NPCNF      : Number of CNFs in AA block                     (Input)
! Cn         : AA Block CI-Coefficients for root selected     (Input)
! lRootSplit : computed root                                  (Input)
! EnFin      : Final Energy for the root selected             (Input)
! DTOC       : Transformation matrix between CSF's and DET's  (input)
! IPRODT     : Prototype determinants                         (input)
! ICONF      : List of configurations                         (input)
! IREFSM     : symmetry of considered CI space                (input)
! Onebod     : one body hamilton matrix in rectangular form   (input)
! ECORE      : Core energy                                    (input)
! NACTOB     : Number of active orbitals                      (input)
! NEL        : total number of active electrons               (input)
! NAEL       : number of alpha active electron                (input)
! NBEL       : number of beta active electron                 (input)
! DIAG       : Hamilton diagonal over CSFs                    (Input)
! TUVX       : Two-electron integrals (MO space)
! IREOTS     : Type => symmetry reordering array
! Ctot       : Vector of all nConf CI-coeff for a single root (Output)

implicit real*8(A-H,O-Z)
#include "spinfo.fh"
#include "WrkSpc.fh"
dimension Ctot(MXPDIM), Cn(NPCSF), IPCSF(MXPDIM), IPCNF(NCONF)
dimension DIAG(*)
dimension DTOC(*), IPRODT(*), ICONF(*)
dimension ONEBOD(*)
dimension TUVX(*), IREOTS(*)

if (NTEST >= 30) then
  write(6,*) ' Input in get_Cm '
  write(6,*) ' ================== '
  write(6,*) ' Total Number of CNFs ',NCONF
  write(6,*) ' Total Number of CSFs ',MXPDIM
  write(6,*) ' CNFs included : '
  call IWRTMA(IPCNF,1,NCONF,1,NCONF)
  write(6,*) ' CSFs included : '
  call IWRTMA(IPCSF,1,MXPDIM,1,MXPDIM)
  write(6,*) ' Number of CNFs in AA block:',NPCNF
  write(6,*) ' Number of CSFs in AA block:',NPCSF
  write(6,*) 'Cn Coefficients'
  call wrtmat(Cn,NPCSF,1,NPCSF,1)
end if

call Fzero(Ctot,MXPDIM)

MXCSFC = 0
do ITYP=1,NTYP
  MXCSFC = max(MXCSFC,NCSFTP(ITYP))
end do

IBblockV = MXPDIM-NPCSF
call getmem('AuxDia','ALLO','REAL',ipAuxD,MXCSFC)
call getmem('AuxGa','ALLO','REAL',ipAuxGa,MXCSFC)
call getmem('AuxGaTi','ALLO','REAL',ipAuxGaTi,MXCSFC)
call getmem('AuxVer','ALLO','REAL',ipAuxV,MXCSFC*NPCSF)
call getmem('AuxBB','ALLO','REAL',ipAuxBB,MXCSFC*IBblockV)
call GETMEM('EXHSCR','MAX','REAL',LW2,MXXWS)
call GETMEM('EXHSCR','ALLO','REAL',LW2,MXXWS)

KACONF = LW2
KLCONF = KACONF+NEL
KLAUXD = KLCONF+NEL
KLFREE = KLAUXD+MXCSFC*MXCSFC

! Shape of matrix:
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
!
! Splitting:
!
!          |
!     AA   |   AB
!  --------|--------
!     BA   |   BB
!          |

call cwtime(C_AlphaLoop1,W_AlphaLoop1)
C_ComputeH_AB = 0.0d0
W_ComputeH_AB = 0.0d0
C_ComputeH_BB = 0.0d0
W_ComputeH_BB = 0.0d0
C_Oper = 0.0d0
W_Oper = 0.0d0

IIAB = 1
do iAlpha=NPCNF+1,NCONF
  !call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
  call cwtime(C_computeH_AB1,W_computeH_AB1)
  if (NTEST >= 30) write(6,*) 'iAlpha = ',iAlpha
  iKACONF = ip_of_iWork_d(Work(KACONF))
  call GETCNF_LUCIA(iWork(iKACONF),IATYP,IPCNF(iAlpha),ICONF,IREFSM,NEL)
  NCSFA = NCSFTP(IATYP)
  if (NTEST >= 30) write(6,*) 'NCSFA = ',NCSFA
  !*********************************************************************
  !                      BB-Block DIAGONAL Elements                    *
  !*********************************************************************
  iKACONF = ip_of_iWork_d(Work(KACONF))
  call CNHCN(iWork(iKACONF),IATYP,iWork(iKACONF),IATYP,Work(KLAUXD),Work(KLFREE),NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX, &
             NTEST,ExFac,IREOTS)
  if (NTEST >= 30) then
    write(6,*) 'Alpha_Alpha elements in BB-block'
    call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
  end if
  do IIA=1,NCSFA
    ILAI = IIA*IIA
    Work(ipAuxD+IIA-1) = Work(KLAUXD+ILAI-1)
    !write(6,*) 'ILAI =',ILAI
    if (NTEST >= 30) then
      write(6,*) 'Work(ipAuxD+IIA-1)',Work(ipAuxD+IIA-1)
    end if
  end do

  IILB = 1
  do Mindex=1,NPCNF ! Loop over AB-Block
    !call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
    if (NTEST >= 30) then
      write(6,*) 'Mindex in AB-Block',Mindex
    end if
    iKLCONF = ip_of_iWork_d(Work(KLCONF))
    call GETCNF_LUCIA(iWork(iKLCONF),ILTYP,IPCNF(Mindex),ICONF,IREFSM,NEL)
    NCSFL = NCSFTP(ILTYP)
    if (NTEST >= 30) write(6,*) 'NCSFL = ',NCSFL
    iKACONF = ip_of_iWork_d(Work(KACONF))
    iKLCONF = ip_of_iWork_d(Work(KLCONF))
    call CNHCN(iWork(iKACONF),IATYP,iWork(iKLCONF),ILTYP,Work(KLAUXD),Work(KLFREE),NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX, &
               NTEST,ExFac,IREOTS)
    if (NTEST >= 30) then
      write(6,*) 'M_Alpha elements'
      call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
    end if
    do IIA=1,NCSFA
      do IIL=1,NCSFL
        IILACT = IILB-1+IIL
        IIAACT = NPCSF*(IIA-1)
        ILAI = (IIL-1)*NCSFA+IIA
        !ILAI = (IIA-1)*MXCSFC+IIL
        ILAOV = IILACT+IIAACT
        Work(ipAuxV+ILAOV-1) = Work(KLAUXD+ILAI-1)
        if (NTEST >= 30) then
          write(6,*) 'ILAI, ILAOV =',ILAI,ILAOV
          write(6,*) 'Work(ipAuxV+ILAOV-1)',Work(ipAuxV+ILAOV-1)
        end if
      end do
    end do
    IILB = IILB+NCSFL
  end do ! End loop over AB-Block
  call cwtime(C_computeH_AB2,W_computeH_AB2)
  C_ComputeH_AB = C_ComputeH_AB+C_computeH_AB2-C_computeH_AB1
  W_ComputeH_AB = W_ComputeH_AB+W_computeH_AB2-W_computeH_AB1

  if (NTEST >= 30) then
    write(6,*) 'AB-Block Vertical Vector'
    call wrtmat(Work(ipAuxV),NPCSF,NCSFA,NPCSF,NCSFA)
  end if

  !*********************************************************************
  ! Compute:                                                           *
  ! 1.  G_Alpha_tilde                                                  *
  ! 2.  G_Alpha                                                        *
  !*********************************************************************
  call cwtime(C_oper1,W_oper1)
  do IIA=1,NCSFA
    Work(ipAuxGaTi+IIA-1) = ddot_(NPCSF,Work(ipAuxV+(IIA-1)*NPCSF),1,Cn,1)
    Work(ipAuxGa+IIA-1) = Work(ipAuxGaTi+IIA-1)/(EnFin-Work(ipAuxD+IIA-1))
    if (NTEST >= 30) then
      write(6,*) 'Work(ipAuxGaTi+IIA-1)',Work(ipAuxGaTi+IIA-1)
      write(6,*) 'Work(ipAuxGa  +IIA-1)',Work(ipAuxGa+IIA-1)
    end if
  end do
  call cwtime(C_oper2,W_oper2)
  C_oper = C_oper+C_oper2-C_oper1
  W_oper = W_oper+W_oper2-W_oper1

  !*********************************************************************
  ! Loop over BB block to get the first order correction to Cm.        *
  ! If not required by the USER (with a keyword) it will be skipped!   *
  !*********************************************************************

  xmaxGaTi = 0.0d0
  do IIA=1,NCSFA
    xmaxGaTi = max(xmaxGaTi,abs(Work(ipAuxGaTi+IIA-1)))
  end do
  if (xmaxGaTi >= 1.0d-12) then ! let's try this trick to make SplitCAS faster!
    call cwtime(C_computeH_BB1,W_computeH_BB1)
    call Fzero(Work(ipAuxBB),MXCSFC*IBblockV)
    IILB = 1
    do Mindex=NPCNF+1,NCONF ! Loop over BB-Block
      !call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
      !if (NTEST >= 30) write(6,*) 'Mindex in BB-Block',Mindex
      iKLCONF = ip_of_iWork_d(Work(KLCONF))
      call GETCNF_LUCIA(iWork(iKLCONF),ILTYP,IPCNF(Mindex),ICONF,IREFSM,NEL)
      NCSFL = NCSFTP(ILTYP)
      !if (NTEST >= 30) write(6,*) 'NCSFL = ',NCSFL
      iKACONF = ip_of_iWork_d(Work(KACONF))
      iKLCONF = ip_of_iWork_d(Work(KLCONF))
      call CNHCN(iWork(iKACONF),IATYP,iWork(iKLCONF),ILTYP,Work(KLAUXD),Work(KLFREE),NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB, &
                 TUVX,NTEST,ExFac,IREOTS)
      !if (NTEST >= 30) then
      !  write(6,*) 'M_Alpha elements in BB-block'
      !  call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
      !end if
      do IIA=1,NCSFA
        do IIL=1,NCSFL
          IILACT = IILB-1+IIL
          !IIAACT = NPCSF*(IIA-1)
          IIAACT = iBblockV*(IIA-1)
          ILAI = (IIL-1)*NCSFA+IIA
          !ILAI = (IIA-1)*MXCSFC+IIL
          ILAOV = IILACT+IIAACT
          Work(ipAuxBB+ILAOV-1) = Work(KLAUXD+ILAI-1)
          if ((iAlpha == Mindex) .and. (IIA == IIL)) then
            !write(6,*) 'iAlpha, IIAB',iAlpha,IIAB
            Work(ipAuxBB+ILAOV-1) = 0.0d0
          end if
          !Write(6,*) 'ILAI, ILAOV =',ILAI,ILAOV
          !if (NTEST >= 30) then
          !  write(6,*) 'Work(ipAuxBB+ILAOV-1)',Work(ipAuxBB+ILAOV-1)
          !end if
        end do
      end do
      IILB = IILB+NCSFL
    end do ! End Loop over BB-Block
    !if (NTEST >= 30) then
    !  write(6,*) 'BB-Block Vertical Vector'
    !  call wrtmat(Work(ipAuxBB),iBblockV,NCSFA,iBblockV,NCSFA)
    !end if
    call cwtime(C_computeH_BB2,W_computeH_BB2)
    C_ComputeH_BB = C_ComputeH_BB+C_computeH_BB2-C_computeH_BB1
    W_ComputeH_BB = W_ComputeH_BB+W_computeH_BB2-W_computeH_BB1

    call cwtime(C_oper1,W_oper1)
    do IIA=1,NCSFA
      call dscal_(iBblockV,Work(ipAuxGa+IIA-1),Work(ipAuxBB+(IIA-1)*iBblockV),1)
      if (NTEST >= 30) then
        write(6,*) 'BB-Block Vertical Vector times Ga'
        call wrtmat(Work(ipAuxBB),iBblockV,NCSFA,iBblockV,NCSFA)
      end if
      call daxpy_(iBblockV,1.0d0,Work(ipAuxBB+(IIA-1)*iBblockV),1,Ctot(NPCSF+1),1)
      if (NTEST >= 30) then
        write(6,*) 'Ctot correction'
        call wrtmat(Ctot,MXPDIM,1,MXPDIM,1)
      end if
    end do
  end if ! End of trick to make SplitCAS faster

  do IIA=1,NCSFA
    Ctot(NPCSF+IIAB+IIA-1) = Ctot(NPCSF+IIAB+IIA-1)+Work(ipAuxGaTi+IIA-1)
    if (NTEST >= 30) then
      write(6,*) 'Ctot '
      call wrtmat(Ctot,MXPDIM,1,MXPDIM,1)
    end if
  end do
  call cwtime(C_oper2,W_oper2)
  C_oper = C_oper+C_oper2-C_oper1
  W_oper = W_oper+W_oper2-W_oper1

  IIAB = IIAB+NCSFA
end do ! End of the loop over iAlpha
if (NTEST >= 30) then
  call cwtime(C_AlphaLoop2,W_AlphaLoop2)
  write(6,*) 'Total time needed to get_Cm in Alpha Loop'
  write(6,*) 'CPU timing : ',C_AlphaLoop2-C_AlphaLoop1
  write(6,*) 'W. timing  : ',W_AlphaLoop2-W_AlphaLoop1

  write(6,*) 'Total time to read H_AB :'
  write(6,*) 'CPU timing : ',C_ComputeH_AB
  write(6,*) 'W. timing  : ',W_ComputeH_AB

  write(6,*) 'Total time to read H_BB :'
  write(6,*) 'CPU timing : ',C_ComputeH_BB
  write(6,*) 'W. timing  : ',W_ComputeH_BB

  write(6,*) 'Total time to calculate (ddot+dscal+daxpy) :'
  write(6,*) 'CPU timing : ',C_Oper
  write(6,*) 'W. timing  : ',W_Oper
end if

call cwtime(C_last1,W_last1)
IIAB = 1
do Mindex=NPCNF+1,NCONF
  !call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
  if (NTEST >= 30) write(6,*) 'Mindex last do loop = ',Mindex
  iKACONF = ip_of_iWork_d(Work(KACONF))
  call GETCNF_LUCIA(iWork(iKACONF),IATYP,IPCNF(Mindex),ICONF,IREFSM,NEL)
  NCSFA = NCSFTP(IATYP)
  if (NTEST >= 30) write(6,*) 'NCSFA = ',NCSFA
  iKACONF = ip_of_iWork_d(Work(KACONF))
  call CNHCN(iWork(iKACONF),IATYP,iWork(iKACONF),IATYP,Work(KLAUXD),Work(KLFREE),NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX, &
             NTEST,ExFac,IREOTS)
  if (NTEST >= 30) then
    write(6,*) 'Alpha_Alpha elements in BB-block'
    call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
  end if
  do IIA=1,NCSFA
    ILAI = IIA*IIA
    !Work(ipAuxD+IIA-1) = Work(KLAUXD+ILAI-1)
    if (NTEST >= 30) then
      write(6,*) 'Work(KLAUXD+ILAI-1)',Work(KLAUXD+ILAI-1)
    end if
    !Work(ipAuxD+IIA-1) = 1.0d0/(EnFin-Work(KLAUXD+ILAI-1))
    Ctot(NPCSF+IIAB+IIA-1) = Ctot(NPCSF+IIAB+IIA-1)/(EnFin-Work(KLAUXD+ILAI-1))
  end do
  IIAB = IIAB+NCSFA
end do
if (NTEST >= 30) then
  call cwtime(C_last2,W_last2)
  write(6,*) 'Time last iteration over Hbb :'
  write(6,*) 'CPU timing : ',C_last2-C_last1
  write(6,*) 'W. timing  : ',W_last2-W_last1
end if

call dcopy_(NPCSF,Cn,1,Ctot,1)
if (NTEST >= 30) then
  write(6,*) 'final Ctot vector'
  call wrtmat(Ctot,MXPDIM,1,MXPDIM,1)
end if
call GETMEM('EXHSCR','FREE','REAL',LW2,MXXWS)
call getmem('AuxBB','FREE','REAL',ipAuxBB,MXCSFC*IBblockV)
call getmem('AuxVer','FREE','REAL',ipAuxV,MXCSFC*NPCSF)
call getmem('AuxGaTi','FREE','REAL',ipAuxGaTi,MXCSFC)
call getmem('AuxGa','FREE','REAL',ipAuxGa,MXCSFC)
call getmem('AuxDia','FREE','REAL',ipAuxD,MXCSFC)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(lrootsPlit)
  call Unused_real_array(DIAG)
end if

end subroutine get_Cm_
