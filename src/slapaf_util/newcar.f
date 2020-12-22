************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine NewCar(Iter,nAtom,nInter,Coor,mTtAtm,Error)
      use Slapaf_Info, only: Cx, qInt, RefGeo, BMx, Shift, Degen,
     &                       AtomLbl, Lbl
      use Slapaf_Parameters, only: Curvilinear, User_Def, BSet, HSet,
     &                             lOld
      Implicit Real*8 (a-h,o-z)
************************************************************************
*                                                                      *
*     Object: To compute the new symm. distinct Cartesian coordinates  *
*             from the suggested shift of the internal coordinates.    *
*                                                                      *
************************************************************************
#include "real.fh"
#include "stdalloc.fh"
#include "weighting.fh"
#include "sbs.fh"
#include "print.fh"
#include "Molcas.fh"
#include "warnings.fh"
      Parameter(NRHS=1)
      Integer, Intent(In):: Iter, nAtom, nInter
      Real*8,  Intent(InOut):: Coor(3,nAtom)
      Integer, Intent(In):: mTtAtm
      Logical, Intent(InOut):: Error
*
      Logical Invar
      Logical:: BSet_Save, HSet_Save, lOld_Save
      Real*8, Allocatable:: DFC(:), dss(:), rInt(:)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt('NewCar: q',' ',qInt,nInter,iter+1)
      Call RecPrt('NewCar: Shift',' ',Shift,nInter,iter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(DFC,3*nAtom,Label='DFC')
      Call mma_allocate(dss,nInter,Label='dss')
      Call mma_allocate(rInt,nInter,Label='rInt')
*                                                                      *
************************************************************************
*                                                                      *
      rInt(:) = qInt(:,iter)
      dss(:) = Shift(:,iter)
      rInt(:) = rInt(:) + dss(:)
*                                                                      *
************************************************************************
*                                                                      *
      Lu=6
      iRout = 33
      iPrint=nPrint(iRout)
      If (iPrint.ge.11) Then
         Write (Lu,*)
         Write (Lu,*) ' *** Transforming internal coordinates'
     &             //' to Cartesian ***'
         Write (Lu,*)
         Write (Lu,*) ' Iter  Internal  Error'
      End If
*
      If (iPrint.ge.99) Then
         Write (Lu,*)
         Write (Lu,*) ' In NewCar: Shifts'
         Write (Lu,*)
         Write (Lu,'(1X,A,2X,F10.4)') (Lbl(iInter),dss(iInter),
     &          iInter=1,nInter)
         Call RecPrt(' In NewCar: qInt',' ',qInt,nInter,
     &               Iter+1)
      End If
*
*     Compute the final internal coordinates, plus sign due to the use
*     of forces and not gradients.
*
      rMax = Zero
      iMax = 0
      jter = 0
      Do i = 1, nInter
         If (Abs(dss(i)) .gt. Abs(rMax)) Then
            rMax = dss(i)
            iMax = i
         End If
      End Do
      If (iPrint.ge.11) Then
         If (iMax.ne.0) Then
            Write (Lu,300) jter,Lbl(iMax),rMax
         Else
            Write (Lu,300) jter,'N/A     ',rMax
         End If
      End If
 300  Format (1X,'Iter:',I5,2X,A,1X,E11.4)
*
      If (iPrint.ge.19) Then
         Write (Lu,*)
         Write (Lu,*)' Internal coordinates of the next macro iteration'
         Write (Lu,*)
         Write (Lu,'(1X,A,2X,F10.4)') (Lbl(iInter),rInt(iInter),
     &          iInter=1,nInter)
      End If
      nPrint_33=nPrint(33)
      nPrint_31=nPrint(31)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the new Cartesian coordinates.
*
      iterMx = 50
*
      Do jter = 1, iterMx
*
*--------Compute the Cartesian shift, solve dq = B^T dx!
*
         M = 3*nAtom
         N = nInter
         Call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,dSS,
     &                  DFC)
         Call mma_deallocate(BMx)
*
         If (iPrint.ge.99)
     &     Call PrList('Symmetry Distinct Nuclear Displacements',
     &                 AtomLbl,nAtom,DFC,3,nAtom)
*
*        Compute the RMS in Cartesian Coordinates.
*
         dx2=Zero
         denom=Zero
         ix = 0
         Do iAtom = 1, nAtom
            Do i = 1, 3
              ix = ix + 1
              dx2 = dx2 + Degen(i,iAtom)*DFC(ix)**2
              denom=denom+Degen(i,iAtom)
            End Do
         End Do
         dx_RMS = Sqrt(dx2/denom)
*
*--------Update the symmetry distinct Cartesian coordinates.
*
         Call DaXpY_(3*nAtom,One,DFC,1,Coor,1)
*
*        Dirty fix of zeros
*
         Do iAtom = 1, nAtom
            If (Cx(1,iAtom,Iter).eq.Zero .and.
     &          Abs(Coor(1,iAtom)).lt.1.0D-13) Coor(1,iAtom)=Zero
            If (Cx(2,iAtom,Iter).eq.Zero .and.
     &          Abs(Coor(2,iAtom)).lt.1.0D-13) Coor(2,iAtom)=Zero
            If (Cx(3,iAtom,Iter).eq.Zero .and.
     &          Abs(Coor(3,iAtom)).lt.1.0D-13) Coor(3,iAtom)=Zero
         End Do
*
         Call dcopy_(3*nAtom,Coor,1,Cx(:,:,Iter+1),1)
         If (iPrint.ge.99)
     &      Call PrList('Symmetry Distinct Nuclear Coordinates / Bohr',
     &                   AtomLbl,nAtom,Coor,3,nAtom)
*
*--------Compute new values q and the Wilson B-matrix for the new
*        geometry with the current new set of Cartesian coordinates.
*
         nWndw=1
         BSet_Save=BSet
         HSet_Save=HSet
         lOld_Save=lOld
         BSet=.False.
         HSet=.False.
         lOld=.False.
         Call BMtrx(nAtom,Coor,iter+1,mTtAtm,nQQ,nWndw)
         BSet=BSet_Save
         HSet=HSet_Save
         lOld=lOld_Save
*
*--------Check if the final structure is reached and get the
*        difference between the present structure and the final.
*
         rOld = rMax
         iMax_Old = iMax
         iMax = 1
         rMax = Zero
         Do i = 1, nInter
            dSS(i) = rInt(i)-qInt(i,Iter+1)
            If (Abs(dSS(i)) .gt. Abs(rMax)) Then
               rMax = dSS(i)
               iMax = i
            End If
         End Do
*
*        Convergence based on the RMS of the Cartesian displacements.
*
         If (dx_RMS.lt.1.0D-6) Go To 998
*
         If (iPrint.ge.99) Then
            Write (Lu,*)
            Write (Lu,*) ' Displacement of internal coordinates'
            Write (Lu,*)
            Write (Lu,'(1X,A,2X,F10.4)') (Lbl(iInter),dss(iInter),
     &             iInter=1,nInter)
         End If
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     On input, Error specifies whether an error should be signalled
*     (.True.) or the calculation aborted (.False.)
*     In the former case, the output value indicates if an error has
*     occurred
*
      If (.NOT.Error) Then
         Call WarningMessage(2,'Error in NewCar')
         Write (Lu,*)
         Write (Lu,*) '***********************************************'
         Write (Lu,*) ' ERROR: No convergence in NewCar !             '
         Write (Lu,*) ' Strong linear dependency among Coordinates.   '
         Write (Lu,*) ' Hint: Try to change the Internal Coordinates. '
         Write (Lu,*) '***********************************************'
         If (.NOT.User_Def) Call RecPrt('NewCar: rInt  ','(10F15.10)',
     &                                  rInt,nInter,1)
         If (.NOT.User_Def) Call RecPrt('NewCar: qInt','(10F15.10)',
     &                                  qInt(:,Iter+1),nInter,1)
         Write (Lu,*)
         Call Quit(_RC_NOT_CONVERGED_)
      End If
*                                                                      *
************************************************************************
*                                                                      *
 998  Continue
*
      If (iPrint.ge.6) Then
         Write (Lu,*)
         Write (Lu,'(A,i2,A)')
     &         ' New Cartesian coordinates were found in',
     &         jter,' Newton-Raphson iterations.'
         Write (Lu,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Finally, just to be safe align the new Cartesian structure with
*     the reference structure (see init2.f)
*
      Invar=(iAnd(iSBS,2**7).eq.0).and.(iAnd(iSBS,2**8).eq.0)
      If (WeightedConstraints.and.Invar)
     &   Call Align(Cx(:,:,iter+1),RefGeo,nAtom)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(rInt)
      Call mma_deallocate(dss)
      Call mma_deallocate(DFC)
*                                                                      *
************************************************************************
*                                                                      *
      Error=.False.
      Return
      End
