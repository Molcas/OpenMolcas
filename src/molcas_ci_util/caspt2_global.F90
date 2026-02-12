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

! Global variables of the CASPT2 module
! TODO: move here all variables from CASPT2 common blocks defined in caspt2_module
module caspt2_global

! UNIT numbers:
! IDCIEX, IDTCEX, LUCIEX, LUDMAT, LUDRA, LUDRATOT, LUH0T, LUHLF1, LUHLF2, LUHLF3, LUINTM, LUONEM, LURHS, LUSBT, LUSOLV
!
! thresholds for printing denominators
! dnmThr, cmpThr, cntThr
! sigma-p regularization, sigma_p_epsilon, sigma_p_exponent

! Some gradient stuff
! iVecL: Solution of the Lambda equation
! idSDMat: offset of overlap derivative; can be defined with 11
! if_SSDM: State-dependent DM is used in Fock or not
! jStLag: The state for which derivatives of the Lagrangian is computed. This is equivalent to jState
!
! unit numbers
! LuPT2, LuGAMMA, LuCMOPT2, LuSTD, LuAPT2
!
! gradients and NAC switches
! do_grad, do_nac, do_csf (CSF term in deriv. coup.)
!
! for removing the weired loop
! iStpGrd, LUGRAD
!
! for IPEA
! do_lindep, if_invar (active invariance), IDSAVGRD, idBoriMat, ConvInvar
!
! whether PT2 energy is invariant wrt rotations among inactive and secondary orbitals
! if_invaria
!
! some derivatives of Lagrangian etc.
! CLag, CLagFull, OLag, OLagFull, SLag, WLag, nCLag, nOLag, nWLag
!
! some correlated density matrices
! DPT2_tot, DPT2C_tot, DPT2_AO_tot, DPT2C_AO_tot, DPT2Canti_tot
!
! Fock-related matrices
! FIMO_all, FIFA_all, FIFASA_all
!
! natural <-> quasi-canonical transformation of frozen orbitals
! TraFro
!
! derivative of the weight factor for XDW-CASPT2
! OMGDER
!
! number of CI vectors per batch in mkfg3.f and derfg3.f
! nbuf1_grad

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: iVecL = 7

integer(kind=iwp) :: idBoriMat(8,13) = 0, IDCIEX, IDSAVGRD = 0, idSDMat(8,13) = 0, IDTCEX, iParRHS, iPrGlb, iRoot1 = 0, &
                     iRoot2 = 0, iStpGrd = 1, jStLag = 0, LuAPT2 = 0, LUCIEX, LuCMOPT2 = 0, LUDMAT, LUDRA, LUDRATOT, LuGAMMA = 0, &
                     LUGRAD = 0, LUH0T(4), LUHLF1, LUHLF2, LUHLF3, LUINTM, LUONEM, LuPT2 = 0, LURHS(8), LUSBT, LUSOLV, LuSTD = 0, &
                     MAXBUF = 1, nbuf1_grad = 0, nCLag = 0, NCMO = 0, NDREF = 0, nOLag = 0, NPREF = 0, nStpGrd = 1, &
                     nTasks_grad = 0, NTAT = 0, NTORB = 0, nWLag = 0, sigma_p_exponent = 2, CompressMPS = 0
real(kind=wp) :: cmpThr, cntThr, ConvInvar = Zero, dnmThr, EMP2, imag_shift = Zero, ipea_shift = Zero, real_shift = Zero, &
                 sigma_p_epsilon = Zero
logical(kind=iwp) :: do_csf = .false., do_grad = .false., do_lindep = .false., do_nac = .false., if_invar = .true., &
                     if_invaria = .true., if_equalW = .true., if_SSDM = .false.
integer(kind=iwp), allocatable :: IDSCT(:), iTasks_grad(:), LISTS(:)
real(kind=wp), allocatable :: CLag(:,:), CLagFull(:,:), CMOPT2(:), DMIX(:,:), DPT2_AO_tot(:), DPT2_tot(:), DPT2C_AO_tot(:), &
                              DPT2C_tot(:), DPT2Canti_tot(:), DREF(:), DWGT(:,:), FAMO(:), FIFA(:), FIFA_all(:), FIFASA_all(:), &
                              FIMO(:), FIMO_all(:), HONE(:), OLag(:), OLagFull(:), OMGDER(:,:), PREF(:), SLag(:,:), TAT(:), &
                              TORB(:), TraFro(:), Weight(:), WLag(:)
real(kind=wp), allocatable, target :: CMO_Internal(:)
real(kind=wp), pointer :: CMO(:)
real(kind=wp), allocatable:: PIQK(:), Buff(:)
integer(kind=iwp), allocatable :: IDXB(:)

public :: CLag, CLagFull, CMO, CMO_Internal, CMOPT2, cmpThr, cntThr, ConvInvar, DMIX, dnmThr, do_csf, do_grad, do_lindep, do_nac, &
          DPT2_AO_tot, DPT2_tot, DPT2C_AO_tot, DPT2C_tot, DPT2Canti_tot, DREF, DWGT, EMP2, FAMO, FIFA, FIFA_all, FIFASA_all, FIMO, &
          FIMO_all, HONE, idBoriMat, IDCIEX, IDSAVGRD, IDSCT, idSDMat, IDTCEX, if_invar, if_invaria, if_equalW, if_SSDM, &
          imag_shift, iParRHS, ipea_shift, iPrGlb, iRoot1, iRoot2, iStpGrd, iTasks_grad, iVecL, jStLag, LISTS, LuAPT2, LUCIEX, &
          LuCMOPT2, LUDMAT, LUDRA, LUDRATOT, LuGAMMA, LUGRAD, LUH0T, LUHLF1, LUHLF2, LUHLF3, LUINTM, LUONEM, LuPT2, LURHS, LUSBT, &
          LUSOLV, LUSTD, MAXBUF, nbuf1_grad, nCLag, NCMO, NDREF, nOLag, NPREF, nStpGrd, NTAT, nTasks_grad, NTORB, nWLag, OLag, &
          OLagFull, OMGDER, PREF, real_shift, sigma_p_epsilon, sigma_p_exponent, SLag, TAT, TORB, TraFro, Weight, WLag, &
          CompressMPS, PIQK, Buff, IDXB

end module caspt2_global
