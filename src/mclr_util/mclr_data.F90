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

module MCLR_Data
! The data of MCLR. The data of MCLR's copy of Lucia (its CI library) are
! in MCLR_Lucia_Data, which only the Lucia routines of mclr_util use.

! MXINKA : Resolution of identity
! MXPORB : Maximum number of orbitals
! Note : MXPNGAS = MXPR4T+6 !!
! Required in order to handle GAS and RAS within /LUCINP/

! Stuff from Pointers.fh
!  ipCI, ipCM, ipMat, ipmatba, ipMatLT, ipMO, n1Dens, n2Dens, na, nacpar, nacpr2, nb, nCMO, nconf1, nDens, nDensC, nmba, nna

! Stuff from genop.fh
! Type of operator in action
!  I12    : 2:Both on and two electron integrals
!  IST    :
!  Square : Integrals square/triangular

! Stuff from machine.fh
!  nrec

! Stuff from csfsd.fh
!  i1, iAllo, iAnders, lconf, lldet

! Stuff from incdia.fh
!  ipdia

! Stuff from sa.fh
!  esterr, Fancy_Preconditioner, irlxroot, isMECIMSPD, isNAC, istate, NACstates, NSSA, override, SA

! Stuff from MCLR_Data.fh
!  ChDisp, dspvec, lDisp, Mxdccc, nhess, SwLbl

! Stuff from MCLR_Data.fh
!  iRefSm  : Reference symmetry
!  MS2     : 2*MS
!  IDC     : Ms combinations
!  PSSIGN  : Ms combination PS factor

! Stuff from Files_mclr.fh
!  FNCSF2SD, FnHlf2, FNHlf3, FnJob, FnMck, FnMOTRA, FnPT2, FnQDAT, FnTemp, FNTrI1, FNTrI2, FNTRi3, FNTRI4, FNTRI5, FnTwo
!  LUCSF2SD, LUHlf2, LUHlf3, LuJob, LuMck, LuMOTRA, LuPt2, LuQDAT, LuTemp, LUTrI1, LUTrI2, LUTRi3, LUTRI4, LUTRI5, LuTwo
!  Nofile

! Stuff from cicisp_mclr.fh
!  niCISp : Number of CI-spaces (1) *
!  iAStFI : Alpha string type for internal CI space *
!  iBStFI : Beta string type for internal CI space (iBZTP)*
!  iActi  : Active/Inactive space  (1)*
!  MnR1IC : Minimum number of electrons in RAS1*
!  MxR1IC : Maximum    -          "      - RAS1*
!  MnR3IC : Minimum    -          "      - RAS3*
!  MxR3IC : Maximum    -          "      - RAS3*
!  nElCI  : Number of electrons per CI space *
!  nAElCI : Number of alpha electrons per CI space *
!  nBElCI : Number of beta electrons per CI space *
!  xispsm : Number of det. for each (symmetry,CI space) **
!  MXSB   : Largest symmetry block **
!  MXSOOB : Largest block        **

! Stuff from crun_mclr.fh
!  iCIStr  : Storing mode for CI-vector
!  NoCSF   : Transform to CSF base
!  iDiag   : PICO
!  Maxi    : Max number of N str. treated simult.(not in use)y
!  MaxK    : Largest number of inner resolution strings simult
!  NoPart  : Partitioning of strings in alpha-alpha + beta-beta loops
!  NACOB   : Active orbiatls
!  NOCOB   : Occupied
!  NTOOB   : Total
!  NORB1   : RAS1
!  NORB2   : RAS2
!  NORB3   : RAS3
!  IREOTS  : Reordering array type     => symmetry
!  ISMFTO  : Symmetry array for type ordered orbitals
!  IBSO    : First orb of given symmetry (symmetry ordered)
!  NTSOB   : Number of active orb of give RAS type and SYM
!  IBTSOB  : Offset for active orb of given RAS type and SYM
!  ITSOB   : Orbitals of given RAS type and sym
!  NOBPTS  : Number of orbitals per subtype and symmetry
!  ITOBS   : "Increaser"
!  ITFSO   : Type array for symmetry ordered orbitals (not activated)

! Stuff from spinfo_mclr.fh
!  MULTSP : Spin multiplicity
!  MS2P   : 2*MS
!  MINOP  : Minum open orbitals     (DONE)
!  MAXOP  : Maximum open orbitals
!  NTYP   : Maxop-MinOp+1           (DONE)
!  NDPCNT : Number of det. / conf type
!  NCPCNT : Number of csf / conf type
!  NCNATS :                         (Done)
!  NDTASM : Combinations
!  NCSASM : CSF
!  NCNASM : i Configurations

! Stuff from CMSLag
!  ResQaaLag2

! Stuff from negpre
!  ERAS, luciv, MXSTATE, ngp, P1, P1INV, SS

! Stuff from PDFT_UTIL
!  Do_Hybrid, PDFT_Ratio, WF_Ratio

! Stuff from Arrays
!  F0SQMO, FAMO, FIMO, Hss
! Various one- and two-particle densities
!  G1t, G2sq, G2t
! MO coefficients
!  CMO, CMO_Inv(:)
!
!  INT1  : 1-electron integrals
!  INT2  : 2-electron integrals
!  PINT1 : Offsets to symmetry blocks
!  KAIN1
!  KINT2
!  KINT2A

! Stuff from Exp
!  H0F, H0S, NewPre, nexp, nexp_max, SBIDT

use Molcas, only: LenIn
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: Mxdccc = 500, MXSTATE = 10


integer(kind=iwp) :: dspvec(mxdccc), ipCI, ipCM(8), ipdia, ipMat(8,8), ipmatba(8,8), ipMatLT(8,8), ipMO(8,8,8), irlxroot, &
                     istate, lDisp(8), luciv, LUHlf2, LUHlf3, LuJob, LuMck, LuMOTRA, LuPt2, LuQDAT, LuTemp, LUTrI1, LUTrI2, &
                     LUTRi3, LUTRI4, LUTRI5, LuTwo, n1Dens, n2Dens, na(8), nacpar, nacpr2, NACstates(2), nb(8), nCMO, nconf1, &
                     nDens, nDensC, nexp = 0, nexp_max = 100, nhess, nmba, nna, nrec, NSSA(2)
real(kind=wp) :: ERAS(MXSTATE), P1(MXSTATE*(MXSTATE+1)/2), P1INV(MXSTATE*(MXSTATE+1)/2), PDFT_Ratio, ResQaaLag2, WF_Ratio
logical(kind=iwp) :: Do_Hybrid, esterr, Fancy_Preconditioner, isMECIMSPD, isNAC, NewPre = .true., ngp, Nofile, override, SA
character(len=LenIn+6) :: ChDisp(Mxdccc*3)
character(len=8) :: FnHlf2, FNHlf3, FnJob, FnMck, FnMOTRA, FnPT2, FnQDAT, FnTemp, FNTrI1, FNTrI2, FNTRi3, FNTRI4, FNTRI5, FnTwo, &
                    SwLbl(Mxdccc)
integer(kind=iwp), allocatable :: H0F(:), SBIDT(:)
real(kind=wp), allocatable :: CMO_Inv(:), F0SQMO(:), FAMO(:), FIMO(:), G1t(:), G2sq(:), G2t(:), H0S(:), Hss(:), INT2(:), SS(:,:)
real(kind=wp), allocatable, target :: CMO(:), INT1(:)

public :: ChDisp, CMO, CMO_Inv, Do_Hybrid, dspvec, ERAS, esterr, F0SQMO, FAMO, Fancy_Preconditioner, FIMO, FnHlf2, FNHlf3, &
          FnJob, FnMck, FnMOTRA, FnPT2, FnQDAT, FnTemp, FNTrI1, FNTrI2, FNTRi3, FNTRI4, FNTRI5, FnTwo, G1t, G2sq, G2t, H0F, H0S, &
          Hss, INT1, INT2, ipCI, ipCM, ipdia, ipMat, ipmatba, ipMatLT, ipMO, irlxroot, isMECIMSPD, isNAC, istate, lDisp, luciv, &
          LUHlf2, LUHlf3, LuJob, LuMck, LuMOTRA, LuPt2, LuQDAT, LuTemp, LUTrI1, LUTrI2, LUTRi3, LUTRI4, LUTRI5, LuTwo, n1Dens, &
          n2Dens, na, nacpar, nacpr2, NACstates, nb, nCMO, nconf1, nDens, nDensC, NewPre, nexp, nexp_max, ngp, nhess, nmba, nna, &
          Nofile, nrec, NSSA, override, P1, P1INV, PDFT_Ratio, ResQaaLag2, SA, SBIDT, SS, SwLbl, WF_Ratio

end module MCLR_Data
