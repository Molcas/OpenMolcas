!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Nike Dattani                                     *
!***********************************************************************

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine level_rdinp(IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,ARV,EPS,NTP,LPPOT,IOMEG1,VLIM1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB, &
                       DSCM,REQ,RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV1,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF)

!***********************************************************************
!  Objective: Read input and construct default input parameters        *
!  written by Nike Dattani in November 2022                            *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: IAN1, IMN1, IAN2, IMN2, CHARGE, NUMPOT, NTP, LPPOT, IOMEG1, IPOTL, PPAR, QPAR, NSR, NLR, IBOB, &
                                  NCMM, IVSR, IDSTT, MMLR(3), NLEV1, AUTO1, LCDC, LXPCT, NJM, JDJR, IWR, LPRWF
real(kind=wp), intent(out) :: RH, RMIN, PRV, ARV, EPS, VLIM1, DSCM, REQ, RREF, RHOAB, CMM(3), PARM(4)
integer(kind=iwp) :: LuIn
logical(kind=iwp) :: skip
character(len=180) :: Line
character(len=4) :: word
character(len=*), parameter :: tabinp(40) = ['IAN1','IMN1','IAN2','IMN2','CHAR','NUMP','RH  ','RMIN','PRV ','ARV ', &
                                             'EPS ','NTP ','LPPO','IOME','VLIM','IPOT','PPAR','QPAR','NSR ','NLR ', &
                                             'IBOB','DSCM','REQ ','RREF','NCMM','IVSR','IDST','RHOA','MMLR','CMM ', &
                                             'PARM','NLEV','AUTO','LCDC','LXPC','NJM ','JDJR','IWR ','LPRW','END ']
integer(kind=iwp), external :: IsFreeUnit
character(len=180), external :: Get_Ln

LuIn = IsFreeUnit(11)
call SpoolInp(LuIn)

! Set default values to input variables

IAN1 = 3                        ! Integer atomic number for atom 1
IMN1 = 6                        ! Integer mass number for atom 1
IAN2 = 3                        ! Integer atomic number for atom 2
IMN2 = 6                        ! Integer mass number for atom 2
CHARGE = 0                      ! Charge of the molecule
NUMPOT = 1

RH = 9.0e-4_wp                  ! Step size (Delta R) for numerical solution of the ODE
RMIN = 0.225_wp                 ! Minimum R for numerical solution of the ODE
PRV = 1                         ! Surkus parameter for the RV (radial variable) for numerical solution of the ODE
ARV = 5.0_wp                    ! Reference distance for the RV (radial variable) for numerical solution of the ODE
EPS = 2.0e-10_wp                ! Epsilon (convergence criterion for numerical solution of the ODE)

NTP = -1                        ! Number of turning points provided (set it to -1 for analytic potentials)
LPPOT = 0                       !
IOMEG1 = 0                      ! Integer Omega quantum number (angular momentum)
VLIM1 = 0                       ! Value of the potential (V) in the R -> infinity limit

IPOTL = 4                       ! Integer potential model switch. Set it to 4 for an MLR model
PPAR = 3                        ! p parameter in the MLR model
QPAR = 3                        ! q parameter in the MLR model
NSR = 3                         ! N_beta (polynomial order, 0 = constant, 1 = linear, etc.) for the short-range side of the potential
NLR = 3                         ! N_beta (polynomial order, 0 = constant, 1 = linear, etc.) for the long-range side of the potential
IBOB = -1

DSCM = 3.337678701485e2_wp      ! D_e (depth of the potential at equilibrium)
REQ = 4.170010583477_wp         ! R_e (equilibrium R value)
RREF = 8.0_wp                   ! R_ref parameter for the Surkus function

NCMM = 3                        ! Number of long-range constants included in u(r) for an MLR model.
IVSR = -2                       !
IDSTT = 1                       !
RHOAB = 0.54_wp                 ! Constant used for damping functions in the MLR model.

MMLR(1) = 6                     ! Inverse-power for the first long-range u(r) term for an MLR model
CMM(1) = 6.719e6_wp             ! Numerator for the first long-range u(r) term for an MLR model
MMLR(2) = 8                     ! Inverse-power for the second long-range u(r) term for an MLR model
CMM(2) = 1.12635e8_wp           ! Numerator for the second long-range u(r) term for an MLR model
MMLR(3) = 10                    ! Inverse-power for the third long-range u(r) term for an MLR model
CMM(3) = 2.78694e9_wp           ! Numerator for the third long-range u(r) term for an MLR model

PARM(1) = -5.156803528943e-1_wp ! Beta_0 for an MLR potential
PARM(2) = -9.585070416286e-2_wp ! Beta_1 for an MLR potential
PARM(3) = 1.170797201140e-1_wp  ! Beta_2 for an MLR potential
PARM(4) = -2.282814434665e-2_wp ! Beta_3 for an MLR potential

NLEV1 = -999
AUTO1 = 1
LCDC = 2
LXPCT = 0
NJM = 0
JDJR = 1
IWR = 3
LPRWF = 0

! Position input file

rewind(LuIn)
call RdNLst(LuIn,'LEVEL')

!ntit1 = 0
skip = .false.

! Read input data from input file

input: do
  if (skip) then
    skip = .false.
  else
    read(LuIn,'(a)') line
    call Upcase(line)
    if (line(1:1) == '*') cycle input
    word = line(1:4)
    if (Word == '') Word = 'END'
  end if

  select case (word)

    ! Read IAN1,IMN1,IAN2,IMN2,CHAR,NUMP
    case (tabinp(1))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IAN1)

    case (tabinp(2))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IMN1)

    case (tabinp(3))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IAN2)

    case (tabinp(4))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IMN2)

    case (tabinp(5))
      Line = Get_Ln(LuIn)
      call Get_I1(1,CHARGE)

    case (tabinp(6))
      Line = Get_Ln(LuIn)
      call Get_I1(1,NUMPOT)

    ! Read RH, RMIN, pRV, aRV, EPS
    case (tabinp(7))
      Line = Get_Ln(LuIn)
      call Get_F1(1,RH)

    case (tabinp(8))
      Line = Get_Ln(LuIn)
      call Get_F1(1,RMIN)

    case (tabinp(9))
      Line = Get_Ln(LuIn)
      call Get_F1(1,PRV)

    case (tabinp(10))
      Line = Get_Ln(LuIn)
      call Get_F1(1,ARV)

    case (tabinp(11))
      Line = Get_Ln(LuIn)
      call Get_F1(1,EPS)
      !write(u6,*) 'EPS=',EPS

    ! Read NTP, LPPOT, IOMEG, VLIM
    case (tabinp(12))
      Line = Get_Ln(LuIn)
      call Get_I1(1,NTP)
      !write(u6,*) 'NTP=',NTP

    case (tabinp(13))
      Line = Get_Ln(LuIn)
      call Get_I1(1,LPPOT)
      !write(u6,*) 'LPPOT=',LPPOT

    case (tabinp(14))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IOMEG1)
      !write(u6,*) 'IOMEG1=',IOMEG1

    case (tabinp(15))
      Line = Get_Ln(LuIn)
      call Get_F1(1,VLIM1)
      !write(u6,*) 'VLIM1=',VLIM1

    ! Read  IPOTL, PPAR, QPAR, NSR, NLR, IBOB
    case (tabinp(16))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IPOTL)
      !write(u6,*) 'IPOTL=',IPOTL

    case (tabinp(17))
      Line = Get_Ln(LuIn)
      call Get_I1(1,PPAR)
      !write(u6,*) 'PPAR=',PPAR

    case (tabinp(18))
      Line = Get_Ln(LuIn)
      call Get_I1(1,QPAR)
      !write(u6,*) 'QPAR=',QPAR

    case (tabinp(19))
      Line = Get_Ln(LuIn)
      call Get_I1(1,NSR)
      !write(u6,*) 'NSR=',NSR

    case (tabinp(20))
      Line = Get_Ln(LuIn)
      call Get_I1(1,NLR)
      !write(u6,*) 'NLR=',NLR

    case (tabinp(21))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IBOB)
      !write(u6,*) 'IBOB=',IBOB

    ! Read DSCM, REQ, Rref
    case (tabinp(22))
      Line = Get_Ln(LuIn)
      call Get_F1(1,DSCM)
      !write(u6,*) 'DSCM=',DSCM

    case (tabinp(23))
      Line = Get_Ln(LuIn)
      call Get_F1(1,REQ)
      !write(u6,*) 'REQ=',REQ

    case (tabinp(24))
      Line = Get_Ln(LuIn)
      call Get_F1(1,RREF)
      !write(u6,*) 'RREF=',RREF

    ! Read NCMM, IVSR, IDSTT, rhoAB
    case (tabinp(25))
      Line = Get_Ln(LuIn)
      call Get_I1(1,NCMM)
      !write(u6,*) 'NCMM=',NCMM

    case (tabinp(26))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IVSR)
      !write(u6,*) 'IVSR=',IVSR

    case (tabinp(27))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IDSTT)
      !write(u6,*) 'IDSTT=',IDSTT

    case (tabinp(28))
      Line = Get_Ln(LuIn)
      call Get_F1(1,RHOAB)
      !write(u6,*) 'RHOAB=',RHOAB
      !write(u6,*) ''

    ! Read MMLR(1), CMM(1), MMLR(2), CMM(2), MMLR(3), CMM(3)
    case (tabinp(29))
      Line = Get_Ln(LuIn)
      call Get_I(1,MMLR,3)
      !write(u6,*) 'MMLR=',MMLR

    case (tabinp(30))
      Line = Get_Ln(LuIn)
      call Get_F(1,CMM,3)
      !write(u6,*) 'CMM=',CMM

    ! Read PARM(1), PARM(2), PARM(3),PARM(4)
    case (tabinp(31))
      Line = Get_Ln(LuIn)
      call Get_F(1,PARM,4)
      !write(u6,*) 'PARM=',PARM

    ! Read NLEV1, AUTO1, LCDC, LXPCT, NJM, JDJR, IWR, LPRWF
    case (tabinp(32))
      Line = Get_Ln(LuIn)
      call Get_I1(1,NLEV1)
      !write(u6,*) 'NLEV1=',NLEV1

    case (tabinp(33))
      Line = Get_Ln(LuIn)
      call Get_I1(1,AUTO1)
      !write(u6,*) 'AUTO1=',AUTO1

    case (tabinp(34))
      Line = Get_Ln(LuIn)
      call Get_I1(1,LCDC)
      !write(u6,*) 'LCDC=',LCDC

    case (tabinp(35))
      Line = Get_Ln(LuIn)
      call Get_I1(1,LXPCT)
      !write(u6,*) 'LXPCT=',LXPCT

    case (tabinp(36))
      Line = Get_Ln(LuIn)
      call Get_I1(1,NJM)
      !write(u6,*) 'NJM=',NJM

    case (tabinp(37))
      Line = Get_Ln(LuIn)
      call Get_I1(1,JDJR)
      !write(u6,*) 'JDJR=',JDJR

    case (tabinp(38))
      Line = Get_Ln(LuIn)
      call Get_I1(1,IWR)
      !write(u6,*) 'IWR=',IWR

    case (tabinp(39))
      Line = Get_Ln(LuIn)
      call Get_I1(1,LPRWF)
      !write(u6,*) 'LPRWF=',LPRWF

    case (tabinp(40))
      exit input

    case default
      write(u6,*)
      write(u6,*) '******************************************'
      write(u6,*) ' LEVEL Error: Input line not recognized.'
      write(u6,*) ' Input line, in upper case:'
      write(u6,'(a)') line
      write(u6,*) ' Extracted keyword: ',word
      write(u6,*) '******************************************'
      call Quit_OnUserError()
  end select
end do input

close(LuIn)

return

end subroutine level_rdinp
