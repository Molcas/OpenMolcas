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

subroutine PCM_Hss(Hess,nHss)

use PCM_arrays
implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "angstr.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "periodic_table.fh"
dimension Hess(nHss)
character*2 Elements(MxAtom*8)

!
!***********************************************************************
!                                                                      *
! Compute the PCM cavity contributions to Hessian in solution

! Retrieve atomic info
call Get_nAtoms_All(nAtoms)
call Allocate_Work(ipCoor,3*nAtoms)
call Get_Coord_All(Work(ipCoor),nAtoms)
call Get_Name_All(Elements)
call GetMem('ANr','Allo','Inte',ipANr,nAtoms)
do i=1,nAtoms
  do j=0,Num_Elem
    if (PTab(j) == Elements(i)) iWork(ipANr+i-1) = j
  end do
end do
call GetMem('Chrg','Allo','Real',ipChrg,nAtoms)
call Get_dArray('Nuclear charge',Work(ipChrg),nAtoms)

! Allocate space for electric field
nComp = 3
call GetMem('EF_n','Allo','Real',ip_EF_n,nComp*nTs)
call GetMem('EF_e','Allo','Real',ip_EF_e,nComp*nTs)

! Allocate space for the matrix derivative, total charges,
! derivatives of charge, two temporary vectors with tesserae
! dimensions, one vector for the derivative of the potential
! on tesserae (with dim. nts*nat3),
! then the PCM contribution to Hessian
nAt3 = nAtoms*3
call GetMem('DerMat','Allo','Real',ip_DerMat,nTs*nTs)
call GetMem('Qtot','Allo','Real',ip_Qtot,nTs)
call GetMem('Qder','Allo','Real',ip_Qder,nAt3*nTs)
call GetMem('Der1','Allo','Real',ip_Der1,nTs)
call GetMem('Der2','Allo','Real',ip_Der2,nTs)
call GetMem('VDer','Allo','Real',ip_VDer,nAt3*nTs)
call GetMem('HssPCM','Allo','Real',ip_HssPCM,nAt3*nAt3)

call PCM_Cav_Hss(Angstr,nAtoms,nAt3,nTs,nS,Eps,iWork(ipANr),Work(ipCoor),Work(ipChrg),Work(ip_EF_n),Work(ip_EF_e),PCMSph,PCMiSph, &
                 PCM_N,PCMTess,PCM_SQ,Work(ip_Qtot),PCMDM,Work(ip_HssPCM),Work(ip_DerMat),dTes,dPnt,dRad,dCntr,Work(ip_QDer), &
                 Work(ip_Der1),Work(ip_Der2),Work(ip_VDer))
call GetMem('HssPCM','Free','Real',ip_HssPCM,nHss)
call GetMem('DerMat','Free','Real',ip_DerMat,nTs*nTs)
call GetMem('Qtot','Free','Real',ip_Qtot,nTs)
call GetMem('Qder','Free','Real',ip_Qder,nAt3*nTs)
call GetMem('Der1','Free','Real',ip_Der1,nTs)
call GetMem('Der2','Free','Real',ip_Der2,nTs)
call GetMem('VDer','Free','Real',ip_VDer,nAt3*nTs)
!                                                                      *
!***********************************************************************
!                                                                      *
return
! Avoid unused argument warnings
if (.false.) call Unused_real_array(Hess)

end subroutine PCM_Hss
