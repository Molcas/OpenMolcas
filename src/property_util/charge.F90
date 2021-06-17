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
! Copyright (C) Valera Veryazov                                        *
!               Luca De Vico                                           *
!***********************************************************************
!  CHARGE
!
!> @brief
!>   Compute and print Mulliken populations
!> @modified_by V. Veryazov
!> @modified_by L. De Vico
!>
!> @details
!> For a given set of natural orbitals compute the
!> Mulliken population analysis. The final print out
!> reports the populations per center and basis
!> function type as well as the gross atomic charges.
!>
!> - \p iCase = ``2`` used for spin independent case (RHF or RASSCF)
!> - \p iCase = ``3`` same as ``2`` but suitable for rasscf *spin* populations.
!> - \p iCase = ``0`` or ``1`` to calculate alpha and beta populations
!>                    (\p iCase = ``0`` only computes alpha contributions, and
!>                    \p iCase = ``1`` calculates beta, and makes a final print out)
!>
!> The Mulliken population is also used to compute the bond
!> order as \f$ \mathit{BO}_{AB} = \sum_{a\in A} \sum_{b\in B} \mathit{DS}_{ab} \mathit{DS}_{ba} \f$.
!>
!> @param[in] NSYM    Number of irreducible representations
!> @param[in] NBAS    Number of basis functions per irred. rep.
!> @param[in] NAME    Center and function type label per basis function
!> @param[in] CMO     Orbital coefficients
!> @param[in] OCCN    Orbital occupations
!> @param[in] SMAT    Overlap matrix
!> @param[in] iCase   Type of run
!> @param[in] FullMlk Boolean for the type of print
!> @param[in] lSave   Boolean for saving on the Runfile
!***********************************************************************

subroutine CHARGE(NSYM,NBAS,NAME,CMO,OCCN,SMAT,iCase,FullMlk,lSave)
! a temporary clone for CHARGE util

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "Molcas.fh"
character*(LENIN8) NAME(*)
dimension NBAS(NSYM), CMO(*), OCCN(*), SMAT(*)
logical FullMlk, lSave, Reduce_Prt
external Reduce_Prt

iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0

if (iCase /= 0) then
  if (iPL >= 2) then
    write(6,*)
    call CollapseOutput(1,'   Molecular charges:')
    write(6,'(3X,A)') '   ------------------'
    write(6,*)
  end if
end if
MXTYP = 0
do iSym=1,nSym
  MxTyp = MxTyp+nBas(iSym)
end do
call Get_iScalar('Unique atoms',nNUC)
call GetMem('QQ','ALLO','REAL',ipQQ,MXTYP*nNuc)
call FZero(Work(ipQQ),MXTYP*nNuc)
call CHARGE_(NSYM,NBAS,NAME,CMO,OCCN,SMAT,iCase,FullMlk,lSave,MXTYP,Work(ipQQ),nNuc)
call GetMem('QQ','FREE','REAL',ipQQ,MXTYP*nNuc)

if (iCase /= 0) then
  if (iPL >= 2) then
    call CollapseOutput(0,'   Molecular charges:')
    write(6,*)
  end if
end if

return

end subroutine CHARGE
