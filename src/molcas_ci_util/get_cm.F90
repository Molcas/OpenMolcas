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

subroutine get_Cm(IPCSF,IPCNF,MXPDIM,NCONF,NPCSF,NPCNF,Cn,EnFin,DTOC,IPRODT,ICONF,IREFSM,ONEBOD,ECORE,NACTOB,NEL,NAEL,NBEL,TUVX, &
                  NTEST,ExFac,IREOTS,FordSplit,Ctot)
!************* Author : GLMJ *****************
!
! Obtain Cm coefficients out of the AA Block
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
! EnFin      : Final Energy for the root selected             (Input)
! DTOC       : Transformation matrix between CSF's and DET's  (Input)
! IPRODT     : Prototype determinants                         (Input)
! ICONF      : List of configurations                         (Input)
! IREFSM     : symmetry of considered CI space                (Input)
! Onebod     : one body hamilton matrix in rectangular form   (Input)
! ECORE      : Core energy                                    (Input)
! NACTOB     : Number of active orbitals                      (Input)
! NEL        : total number of active electrons               (Input)
! NAEL       : number of alpha active electron                (Input)
! NBEL       : number of beta active electron                 (Input)
! TUVX       : Two-electron integrals (MO space)
! IREOTS     : Type => symmetry reordering array
! Ctot       : Vector of all nConf CI-coeff for a single root (Output)

use spinfo, only: NCSFTP, NTYP
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MXPDIM, IPCSF(MXPDIM), NCONF, IPCNF(NCONF), NPCSF, NPCNF, IPRODT(*), ICONF(*), IREFSM, NACTOB, &
                                 NEL, NAEL, NBEL, IREOTS(NACTOB)
real(kind=wp), intent(in) :: Cn(NPCSF), EnFin, DTOC(*), ONEBOD(NACTOB,NACTOB), ECORE, TUVX(*), ExFac
integer(kind=iwp), intent(inout) :: NTEST
logical(kind=iwp), intent(in) :: FordSplit
real(kind=wp), intent(out) :: Ctot(MXPDIM)
integer(kind=iwp) :: iAlpha, IATYP, IIA, IIAB, IIL, IILACT, IILB, ILAI, ILTYP, ITYP, Mindex, MXCSFC, MXXWS, NCSFA, NCSFL
real(kind=wp) :: C_AlphaLoop1, C_AlphaLoop2, C_ComputeH_AB, C_computeH_AB1, C_computeH_AB2, C_Oper, C_oper1, C_oper2, &
                 W_AlphaLoop1, W_AlphaLoop2, W_ComputeH_AB, W_computeH_AB1, W_computeH_AB2, W_Oper, W_oper1, W_oper2
integer(kind=iwp), allocatable :: ICNL(:), ICNR(:)
real(kind=wp), allocatable :: AuxD(:), AuxGa(:), AuxGaTi(:), AuxV(:,:), CNHCNM(:), Scr(:)
real(kind=wp), external :: ddot_

if (NTEST >= 30) then
  write(u6,*) ' Input in get_Cm'
  write(u6,*) ' =================='
  write(u6,*) ' Total Number of CNFs ',NCONF
  write(u6,*) ' Total Number of CSFs ',MXPDIM
  write(u6,*) ' CNFs included :'
  call IWRTMA(IPCNF,1,NCONF,1,NCONF)
  write(u6,*) ' CSFs included :'
  call IWRTMA(IPCSF,1,MXPDIM,1,MXPDIM)
  write(u6,*) ' Number of CNFs in AA block:',NPCNF
  write(u6,*) ' Number of CSFs in AA block:',NPCSF
  write(u6,*) 'Cn Coefficients'
  call wrtmat(Cn,NPCSF,1,NPCSF,1)
end if
!***********************************************************************
! get Cm coefficients corrected to the firts order                     *
! according to Lowdin's equations                                      *
!***********************************************************************
if (FOrdSplit) then
  call get_Cm_(IPCSF,IPCNF,MXPDIM,NCONF,NPCSF,NPCNF,Cn,EnFin,DTOC,IPRODT,ICONF,IREFSM,ONEBOD,ECORE,NACTOB,NEL,NAEL,NBEL,TUVX, &
               NTEST,ExFac,IREOTS,Ctot)
  return
end if

!***********************************************************************
! get Cm coefficients corrected to the zeroth order                    *
! according to Lowdin's equations                                      *
!***********************************************************************
Ctot(:) = Zero

MXCSFC = 0
do ITYP=1,NTYP
  MXCSFC = max(MXCSFC,NCSFTP(ITYP))
end do

call mma_allocate(AuxD,MXCSFC,label='AuxDia')
call mma_allocate(AuxGa,MXCSFC,label='AuxGa')
call mma_allocate(AuxGaTi,MXCSFC,label='AuxGaTi')
call mma_allocate(AuxV,NPCSF,MXCSFC,label='AuxVer')

call mma_allocate(ICNL,NEL,Label='ICNL')
call mma_allocate(ICNR,NEL,Label='ICNR')
call mma_allocate(CNHCNM,MXCSFC**2,Label='CNHCNM')

call mma_maxDBLE(MXXWS)
call mma_allocate(Scr,MXXWS,label='EXHSCR')

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
C_ComputeH_AB = Zero
W_ComputeH_AB = Zero
C_Oper = Zero
W_Oper = Zero

IIAB = 1
do iAlpha=NPCNF+1,NCONF
  !CNHCNM(:) = Zero
  call cwtime(C_computeH_AB1,W_computeH_AB1)
  if (NTEST >= 30) write(u6,*) 'iAlpha = ',iAlpha
  call GETCNF(ICNL,IATYP,IPCNF(iAlpha),ICONF,IREFSM,NEL)
  NCSFA = NCSFTP(IATYP)
  if (NTEST >= 30) write(u6,*) 'NCSFA = ',NCSFA
  !*********************************************************************
  !                      BB-Block DIAGONAL Elements                    *
  !*********************************************************************
  call CNHCN(ICNL,IATYP,ICNR,IATYP,CNHCNM,SCR,NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,NTEST,ExFac,IREOTS)
  if (NTEST >= 30) then
    write(u6,*) 'Alpha_Alpha elements in BB-block'
    call wrtmat(CNHCNM,MXCSFC,MXCSFC,MXCSFC,MXCSFC)
  end if
  do IIA=1,NCSFA
    ILAI = IIA*IIA
    AuxD(IIA) = CNHCNM(ILAI)
    !write(u6,*) 'ILAI =',ILAI
    if (NTEST >= 30) write(u6,*) 'AuxD(IIA)',AuxD(IIA)
  end do

  IILB = 1
  do Mindex=1,NPCNF ! Loop over AB-Block
    !CNHCNM(:) = Zero
    if (NTEST >= 30) then
      write(u6,*) 'Mindex in AB-Block',Mindex
    end if
    call GETCNF(ICNR,ILTYP,IPCNF(Mindex),ICONF,IREFSM,NEL)
    NCSFL = NCSFTP(ILTYP)
    if (NTEST >= 30) write(u6,*) 'NCSFL = ',NCSFL
    call CNHCN(ICNL,IATYP,ICNR,ILTYP,CNHCNM,SCR,NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,NTEST,ExFac,IREOTS)
    if (NTEST >= 30) then
      write(u6,*) 'M_Alpha elements'
      call wrtmat(CNHCNM,MXCSFC,MXCSFC,MXCSFC,MXCSFC)
    end if
    do IIA=1,NCSFA
      do IIL=1,NCSFL
        IILACT = IILB-1+IIL
        ILAI = (IIL-1)*NCSFA+IIA
        !ILAI = (IIA-1)*MXCSFC+IIL
        AuxV(IILACT,IIA) = CNHCNM(ILAI)
        if (NTEST >= 30) then
          write(u6,*) 'ILAI, IILACT, IIA =',ILAI,IILACT,IIA
          write(u6,*) 'AuxV(IILACT,IIA)',AuxV(IILACT,IIA)
        end if
      end do
    end do
    IILB = IILB+NCSFL
  end do ! End loop over AB-Block
  call cwtime(C_computeH_AB2,W_computeH_AB2)
  C_ComputeH_AB = C_ComputeH_AB+C_computeH_AB2-C_computeH_AB1
  W_ComputeH_AB = W_ComputeH_AB+W_computeH_AB2-W_computeH_AB1

  if (NTEST >= 30) then
    write(u6,*) 'AB-Block Vertical Vector'
    call wrtmat(AuxV,NPCSF,NCSFA,NPCSF,NCSFA)
  end if

  !*********************************************************************
  ! Compute:                                                           *
  ! 1.  G_Alpha_tilde                                                  *
  ! 2.  G_Alpha                                                        *
  !*********************************************************************
  call cwtime(C_oper1,W_oper1)
  do IIA=1,NCSFA
    AuxGaTi(IIA) = ddot_(NPCSF,AuxV(:,IIA),1,Cn,1)
    AuxGa(IIA) = AuxGaTi(IIA)/(EnFin-AuxD(IIA))
    if (NTEST >= 30) then
      write(u6,*) 'AuxGaTi(IIA)',AuxGaTi(IIA)
      write(u6,*) 'AuxGa(IIA)  ',AuxGa(IIA)
    end if
  end do
  call cwtime(C_oper2,W_oper2)
  C_oper = C_oper+C_oper2-C_oper1
  W_oper = W_oper+W_oper2-W_oper1

  do IIA=1,NCSFA
    Ctot(NPCSF+IIAB+IIA-1) = Ctot(NPCSF+IIAB+IIA-1)+AuxGa(IIA)
    if (NTEST >= 30) then
      write(u6,*) 'Ctot'
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
  write(u6,*) 'Total time needed to get_Cm in Alpha Loop'
  write(u6,*) 'CPU timing : ',C_AlphaLoop2-C_AlphaLoop1
  write(u6,*) 'W. timing  : ',W_AlphaLoop2-W_AlphaLoop1

  write(u6,*) 'Total time to read H_AB :'
  write(u6,*) 'CPU timing : ',C_ComputeH_AB
  write(u6,*) 'W. timing  : ',W_ComputeH_AB

  write(u6,*) 'Total time to calculate (ddot+dscal+daxpy) :'
  write(u6,*) 'CPU timing : ',C_Oper
  write(u6,*) 'W. timing  : ',W_Oper
end if
Ctot(1:NPCSF) = Cn(:)
if (NTEST >= 30) then
  write(u6,*) 'final Ctot vector'
  call wrtmat(Ctot,MXPDIM,1,MXPDIM,1)
end if
call mma_deallocate(AuxD)
call mma_deallocate(AuxGa)
call mma_deallocate(AuxGaTi)
call mma_deallocate(AuxV)
call mma_deallocate(ICNL)
call mma_deallocate(ICNR)
call mma_deallocate(CNHCNM)
call mma_deallocate(Scr)

return

end subroutine get_Cm
