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

subroutine read_input(IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,pRV,aRV,EPS,NTP,LPPOT,IOMEG,VLIM,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,Rref,NCMM,IVSR,TDSTT,rhoAB,MMLR,CMM,PARM,NLEV,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF)

!***********************************************************************
!  Objective: Read input and construct default input parameters        *
!  written by Nike Dattani in November 2022                            *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) ::  IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,pRV,NTP,LPPOT,IOMEG,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,NCMM,IVSR,TDSTT,MMLR,NLEV,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF
real(kind=wp), intent(out) :: RH,RMIN,EPS,VLIM,DSCM,REQ,Rref,rhoAB,CMM,PARM
real(kind=wp), intent(out) :: Rout(npoint+4), PotR(npoint+4)
character(len=80) :: Title1(10), Title2(10)
integer(kind=iwp), parameter :: ntab = 40
character(len=*), parameter :: tabinp(ntab) = ['IAN1','IMN1','IAN2','IMN2','CHAR','NUMP','RH  ','RMIN','pRV ','aRV ', &
                                               'EPS ','NTP ','LPPO','IOME','VLIM','IPOT','PPAR','QPAR','NSR ','NLR ', &
                                               'IBOB','DSCM','REQ ','Rref','NCMM','IVSR','TDST','rhoA','MMLR','CMM ', &
                                               'PARM','NLEV','AUTO','LCDC','LXPC','NJM ','JDJR','IWR ','LPRW','END ']
character(len=180), external :: Get_Ln, Get_Ln_EOF

LuIn = IsFreeUnit(11)
call SpoolInp(LuIn)

! Set default values to input variables

IAN1 = 3        ! Integer atomic number for atom 1
IMN1 = 6        ! Integer mass number for atom 1
IAN2 = 3        ! Integer atomic number for atom 2
IMN2 = 6        ! Integer mass number for atom 2
CHARGE = 0      ! Charge of the molecule
NUMPOT = 1      ! 

RH = 0.0005     ! Step size (Delta R) for numerical solution of the ODE
RMIN = 0.125    ! Minimum R for numerical solution of the ODE
pRV = 1         ! Surkus parameter for the RV (radial variable) for numerical solution of the ODE
aRV = 5.0       ! Reference distance for the RV (radial variable) for numerical solution of the ODE
EPS = 2.d-10    ! Epsilon (convergence criterion for numerical solution of the ODE)

NTP = -1        ! Number of turning points provided (set it to -1 for analytic potentials)
LPPOT = 0       ! 
IOMEG = 0       ! Integer Omega quantum number (angular momentum)
VLIM = 0        ! Value of the potential (V) in the R -> infinity limit

IPOTL = 4       ! Integer potential model switch. Set it to 4 for an MLR model
PPAR = 3        ! p parameter in the MLR model
QPAR = 3        ! q parameter in the MLR model
NSR = 3         ! N_beta (polynomial order, 0 for a constant, 1 for linear, etc.) for the short-range side of the potential
NLR = 3         ! N_beta (polynomial order, 0 for a constant, 1 for linear, etc.) for the long-range side of the potential
IBOB = -1       ! 

DSCM = 3.337678701485D+02 ! D_e (depth of the potential at equilibrium)
REQ = 4.170010583477D+00  ! R_e (equilibrium R value)
Rref = 8.0d0              ! R_ref parameter for the Surkus function

NCMM = 3        ! Number of long-range constants included in u(r) for an MLR model.
IVSR  = -2      !
TDSTT = 1       ! 
rhoAB = 0.54    ! Constant used for damping functions in the MLR model.

MMLR(1) = 6          ! Inverse-power for the first long-range u(r) term for an MLR model
CMM(1) = 6.719d+06   ! Numerator for the first long-range u(r) term for an MLR model
MMLR(2) = 8          ! Inverse-power for the second long-range u(r) term for an MLR model
CMM(2) = 1.12635d+08 ! Numerator for the second long-range u(r) term for an MLR model
MMLR(3) = 10         ! Inverse-power for the third long-range u(r) term for an MLR model  
CMM(3) = 2.78694d+09 ! Numerator for the third long-range u(r) term for an MLR model

PARM(1) = -5.156803528943D-01 ! Beta_0 for an MLR potential
PARM(2) = -9.585070416286D-02 ! Beta_1 for an MLR potential
PARM(3) =  1.170797201140D-01 ! Beta_2 for an MLR potential
PARM(4) = -2.282814434665D-02 ! Beta_3 for an MLR potential

NLEV = -999     !
AUTO1 = 1       !
LCDC = 2        !
LXPCT = 0       !
NJM = 0         !
JDJR = 1        !
IWR = 3         !
LPRWF = 0       !

! Position input file

rewind(LuIn)
call RdNLst(LuIn,'LEVEL')

! Read input data from input file

    case (tabinp(1))
      ! Read atomic information
      ! IAN1 and IAN2 are the atomic numbers for atoms 1 and 2,
      ! IMN1 and IMN2 are the corresponding mass numbers
      Line = Get_Ln(LuIn)
      call Get_I1(1,IAN1)

    case (tabinp(7))
      ! Read the step size RH
      Line = Get_Ln(LuIn)
      call Get_F1(1,RH)

    case (tabinp(40))
      exit input

return

end subroutine read_input
