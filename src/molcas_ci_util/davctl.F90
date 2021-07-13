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
! Copyright (C) 1989, Bjorn O. Roos                                    *
!               1989, Per Ake Malmqvist                                *
!               1991, Jeppe Olsen                                      *
!               1991,1996, Markus P. Fuelscher                         *
!               2000, Thorstein Thorsteinsson                          *
!***********************************************************************

subroutine DavCtl(LW1,TUVX,IFINAL)
!***********************************************************************
!                                                                      *
!     CI control section                                               *
!                                                                      *
!     calling arguments:                                               *
!     LW1     : active Fock matrix                                     *
!               array of real*8                                        *
!     TUVX    : array of real*8                                        *
!               two-electron integrals (tu!vx)                         *
!     IFINAL  : integer                                                *
!               termination flag                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     B.O. Roos and P.Aa. Malmqvist                                    *
!     University of Lund, Sweden, 1989                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     - updated to use determinant based CI-procedures                 *
!       J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
!     - updated for MOLCAS version 3                                   *
!       J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
!     - updated for integral direct and reaction field calculations    *
!       M.P. Fuelscher, University of Lund, Sweden, 1996               *
!     - various modifications                                          *
!       T. Thorsteinsson, University of Lund, Sweden, 2000             *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
real*8 LW1
dimension LW1(*), TUVX(*)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "ciinfo.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
#include "lucia_ini.fh"

!-----------------------------------------------------------------------
!    INITIALIZE THE DAVIDSON DIAGONALIZATION
!-----------------------------------------------------------------------

lRoots = lRoots+hroots
call Ini_David(lRoots,nConf,nDet,nSel,n_keep,nAc,LuDavid)

IPRLEV = IPRLOC(3)

!-----------------------------------------------------------------------
!     COMPUTE THE DIAGONAL ELEMENTS OF THE HAMILTONIAN
!-----------------------------------------------------------------------

! LW4: TEMPORARY CI VECTOR IN CSF BASIS

if (IPRLEV >= 20) write(6,1100) 'INI_DAVID'
call GETMEM('CIVEC','ALLO','REAL',LW4,NCONF)
if (IPRLEV >= 20) write(6,1100) 'CIDIA',LW4
if (NAC > 0) call CIDIA_CI_UTIL(NAC,NCONF,STSYM,WORK(LW4),LW1,TUVX,LUDAVID)

!-----------------------------------------------------------------------
!    OBTAIN STARTING VECTORS
!-----------------------------------------------------------------------

mSel = nSel
if (NAC /= 0) then
  call GetMem('iSel','Allo','Integer',lSel,mSel)
  call GetMem('ExplE','Allo','Real',lExplE,mSel)
  call GetMem('ExplV','Allo','Real',lExplV,mSel*mSel)
else
  lSel = 1
  lExplE = 1
  lExplV = 1
end if
nMaxSel = nConf
if (N_ELIMINATED_GAS_MOLCAS > 0) nmaxSel = nCSF_HEXS

if ((IPRLEV >= 20) .and. (NAC /= 0)) write(6,1100) 'CSTART',LW4,lSel,lExplE,lExplV
call CStart_CI_Util(WORK(LW4),LW1,TUVX,iWork(lSel),Work(lExplE),Work(lExplV),nMaxSel,IFINAL)

call GETMEM('CIVEC','FREE','REAL',LW4,NCONF)

!-----------------------------------------------------------------------
!    DIAGONALIZATION SECTION
!-----------------------------------------------------------------------

! PAM Jun 2006: Gradual lowering of threshold in first 3 iterations
!if (Iter > 1) Threshold = THFACT*abs(CONV(4,ITER-1))
! Energy threshold for CI convergence criterion:
if (Iter == 1) then
  Threshold = THREN
else if ((ITER > 1) .and. (ITER <= 3)) then
  ThrRule = THFACT*abs(CONV(4,ITER-1))
  Threshold = (dble(4-ITER)*THREN+dble(ITER)*ThrRule)*0.25d0
else
  Threshold = THFACT*abs(CONV(4,ITER-1))
end if
! End of new rule, PAM Jun 2006
Threshold = max(Threshold,1.0D-9)
if (NAC == 0) then
  ESize = abs(EMY)
else
  ESize = abs(Work(lExplE))
end if
Threshold = max(Threshold,ESize*1.0D-14)

! LW5: CONVERGENCE PARAMETERS
call GetMem('CI_conv','Allo','Real',LW5,2*lRoots*MAXJT)
if ((IPRLEV >= 20) .and. (NAC /= 0)) write(6,1100) 'DAVID',LW5,lSel,lExplE,lExplV
ITERCI = 1
if (NAC == 0) then
  ENER(1,ITER) = EMY
else
  if ((nSel == nConf) .or. ((N_ELIMINATED_GAS_MOLCAS > 0) .and. (nSel == nCSF_HEXS))) then
    do jRoot=1,lRoots-hRoots
      ENER(jRoot,ITER) = Work(lExplE+jRoot-1)
    end do
  else
    ! PAM Jun 2006: Limit the number of CI iterations in the
    ! first few macroiterations:
    if (KTIGHT == 0) ItLimit = min(12*ITER,MAXJT)
    if (KTIGHT == 1) ItLimit = MAXJT
    ! PAM Oct 2006: Full precision if this is final CI.
    if ((ICIONLY == 1) .or. (IFINAL == 2)) then
      Threshold = max(1.0D-9,ESize*1.0D-14)
      ITLIMIT = MAXJT
    end if
    ! PAM Feb 2009: New code in david5.
    !call David5(nAc,stSym,nDet,MAXJT,ITERCI,
    !call David5(nAc,stSym,nDet,ItLimit,ITERCI,Work(LW5),Threshold,LW1,TUVX,iWork(lSel),Work(lExplE),Work(lExplV))

    call David5(nDet,ItLimit,IterCI,Work(LW5),Threshold,iWork(lSel),Work(lExplE),Work(lExplV),LW1,TUVX)

    do jRoot=1,lRoots-hRoots
      ENER(jRoot,ITER) = Work(LW5+2*(jRoot-1)+2*(ITERCI-1)*lRoots)
    end do
  end if
end if
call GetMem('CI_conv','Free','Real',LW5,2*lRoots*MAXJT)
if (NAC /= 0) then
  call GetMem('ExplV','Free','Real',lExplV,mSel*mSel)
  call GetMem('ExplE','Free','Real',lExplE,mSel)
  call GetMem('iSel','Free','Integer',lSel,mSel)
end if
nSel = mSel
lRoots = lRoots-hroots

!-----------------------------------------------------------------------
!    CLEANUP AFTER THE DAVIDSON DIAGONALIZATION
!-----------------------------------------------------------------------

! LW4: TEMPORARY CI VECTOR IN CSF BASIS

call GETMEM('CIVEC','ALLO','REAL',LW4,NCONF)
if (IPRLEV >= 20) write(6,1100) 'TERM_DAVID',LW4
iDisk = IADR15(4)
call Term_David(ICICH,ITERCI,lRoots,nConf,Work(LW4),JOBIPH,LuDavid,iDisk)
call GETMEM('CIVEC','FREE','REAL',LW4,NCONF)

return

1100 format(1X,/,1X,'WORK SPACE VARIABLES IN SUBR. CICTL: ',/,1X,'SUBSECTION: ',A,/,(1X,12I10,/))

end subroutine DavCtl
