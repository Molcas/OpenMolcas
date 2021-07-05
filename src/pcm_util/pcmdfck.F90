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

subroutine PCMDFck(nFck,PCMFck)

use PCM_arrays
implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "Molcas.fh"
#include "rctfld.fh"
#include "periodic_table.fh"
dimension PCMFck(nFck,*)
character*2 Elements(MxAtom*8)
logical DoPot, DoFld

!
!***********************************************************************
!
! Driver for the computation of PCM contributions to Fock matrix
! derivatives

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

! Allocate space for total charges
call GetMem('Qtot','Allo','Real',ip_Qtot,nTs)

! Allocate space for the PCM matrix derivative
call GetMem('DerMat','Allo','Real',ip_DerMat,nTs*nTs)

! Allocate space for the potential on tesserae
call GetMem('V','Allo','Real',ip_V,nTs)

! Allocate space for the uncontracted potential on tesserae
call GetMem('VMN','Allo','Real',ip_VMN,nTs*nFck)

! Allocate space for the potential derivatives
nAt3 = nAtoms*3
call GetMem('VDer','Allo','Real',ip_VDer,nAt3*nTs)

! Allocate space for the uncontracted derivatives of the
! potential on tesserae
call GetMem('VDerMN','Allo','Real',ip_VDerMN,nAt3*nTs*nFck)

! Allocate space for electric field
nComp = 3
call GetMem('EF_n','Allo','Real',ip_EF_n,nComp*nTs)
call GetMem('EF_e','Allo','Real',ip_EF_e,nComp*nTs)

! Allocate two scratch vectors
call GetMem('Temp1','Allo','Real',ip_Temp1,nTs)
call GetMem('Temp2','Allo','Real',ip_Temp2,nTs)

call FZero(PCMFck,nAt3*nFck)

! Compute the potential and the electric field on tesserae
DoPot = .true.
DoFld = .true.
call V_EF_PCM(nAtoms,nTs,DoPot,DoFld,Work(ipCoor),PCMTess,Work(ip_V),Work(ip_EF_n),Work(ip_EF_e))

! Compute the derivatives of the total potential on tesserae
call VDer_PCM(nAtoms,nTs,nS,Work(ipCoor),Work(ipChrg),Work(ip_EF_n),Work(ip_EF_e),PCMTess,PCMiSph,dTes,dPnt,dRad,dCntr, &
              Work(ip_VDer))

! Actually compute the PCM correction
call PCM_Der_Fock(nFck,nAtoms,nTs,nS,Eps,PCMSph,PCMiSph,PCM_N,PCMTess,PCM_SQ,Work(ip_Qtot),PCMDM,Work(ip_DerMat),dTes,dPnt,dCntr, &
                  Work(ip_V),Work(ip_VMN),Work(ip_VDer),Work(ip_VDerMN),Work(ip_Temp1),Work(ip_Temp2),PCMFck)

! Free the space
call GetMem('Qtot','Free','Real',ip_Qtot,nTs)
call GetMem('V','Free','Real',ip_V,nTs)
call GetMem('VMN','Free','Real',ip_VMN,nTs*nFck)
call GetMem('VDer','Free','Real',ip_VDer,nAt3*nTs)
call GetMem('VDerMN','Free','Real',ip_VDerMN,nAt3*nTs*nFck)
call GetMem('EF_n','Free','Real',ip_EF_n,nComp*nTs)
call GetMem('EF_e','Free','Real',ip_EF_e,nComp*nTs)
call GetMem('Temp1','Free','Real',ip_Temp1,nTs)
call GetMem('Temp2','Free','Real',ip_Temp2,nTs)

return

end subroutine PCMDFck
!====
subroutine PCM_Der_Fock(nFck,nAt,nTs,nS,Eps,Sphere,ISphe,nOrd,Tessera,Q,Qtot,DM,DerDM,DerTes,DerPunt,DerCentr,V,VMN,VDer,VDerMN, &
                        Temp1,Temp2,PCMFck)

implicit real*8(a-h,o-z)
#include "real.fh"
dimension Sphere(4,*), ISphe(*), nOrd(*), Temp1(*), Temp2(*)
dimension Tessera(4,*), Q(2,*), Qtot(*), DM(nTs,*), DerDM(nTs,*)
dimension V(*), VMN(nTs,*)
dimension VDer(nTs,*), VDerMN(nTs,nFck,*), PCMFck(nFck,*)
dimension DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3)
dimension DerCentr(nS,nAt,3,3)

FPI = Four*PI
Diag = -1.0694d0*sqrt(FPI)/Two

do iTs=1,nTs
  Qtot(iTs) = Q(1,iTs)+Q(2,iTs)
end do

! Loop over the degrees of freedom
do iAt=1,nAt
  do iC=1,3
    Index = 3*(iAt-1)+iC

    ! Derivative of the PCM matrix for the conductor-like case
    call DMat_CPCM(iAt,iC,Eps,nTs,nS,nAt,Diag,Tessera,DerDM,DerTes,DerPunt,DerCentr,iSphe)

    ! Solvation charges (weights) times the derivative of the PCM matrix
    call PrMatVec(.true.,.true.,DerDM,-1.d0,nTs,nTs,QTot,Temp1)

    ! The previous vector times the inverted PCM matrix
    call PrMatVec(.true.,.true.,DM,1.d0,nTs,nTs,Temp1,Temp2)

    ! Derivative of the potential times the inverted PCM matrix
    call PrMatVec(.true.,.true.,DM,1.d0,nTs,nTs,VDer(1,Index),Temp1)

    ! Loop over the Fock elements
    do iFck=1,nFck

      ! First contribution: charges times the derivative of the
      ! (uncontracted) electronic potential: Sum_i q_i V_mn^(i,x)
      Sum = Zero
      do iTs=1,nTs
        Sum = Sum+Qtot(iTs)*VDerMN(iTs,iFck,Index)
      end do
      PCMFck(iFck,Index) = PCMFck(iFck,Index)+Sum

      ! Second contribution: derivative of the potential times the
      ! (symmetrized) PCM matrix times the uncontracted potential:
      ! Sum_ij V_i^x Q_ij V_mn^j
      Sum = Zero
      do iTs=1,nTs
        Sum = Sum+Temp1(iTs)*VMN(iTs,iFck)
      end do
      PCMFck(iFck,Index) = PCMFck(iFck,Index)+Sum

      ! Third contribution: potential times the derivative of
      ! the inverted PCM matrix times the uncontracted potential:
      ! Sum_ij V_i Q_ij^x V_mn^j = Sum_ij V_i [-S^-1 S^x S^-1]_ij V_mn^j =
      ! Sum_ij q_i [S^x S^-1]_ij V_mn^j
      Sum = Zero
      do iTs=1,nTs
        do jTs=1,nTs
          Sum = Sum+Temp2(iTs)*VMN(jTs,iFck)
        end do
      end do
      PCMFck(iFck,Index) = PCMFck(iFck,Index)+Sum
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Sphere)
  call Unused_integer_array(nOrd)
  call Unused_real_array(V)
end if

end subroutine PCM_Der_Fock
