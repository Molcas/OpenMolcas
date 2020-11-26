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
      Subroutine Int2Car(dSS,rInt,nInter,ip_qInt,Coor,nAtom,nBVct,
     &                  ipBMx,dMass,nLines,DFC,
     &                  nDim,Lbl,Name,iSym,
     &                  Smmtrc,Degen,iter,
     &                  ip_dqInt,mTtAtm,iANr,iOptH,
     &                  User_Def,nStab,jStab,Curvilinear,
     &                  Numerical,DDV_Schlegel,HWRS,
     &                  Analytic_Hessian,iOptC,PrQ,mxdc,
     &                  iCoSet,rHidden,Error,ipRef,Redundant,nqInt,
     &                  MaxItr,iRef)
      use Slapaf_Info, only: Cx, Gx
      Implicit Real*8 (a-h,o-z)
************************************************************************
*                                                                      *
*     Object: To compute the new symm. distinct Cartesian coordinates  *
*             from the suggested shift of the internal coordinates.    *
*                                                                      *
************************************************************************
#include "real.fh"
#include "WrkSpc.fh"
#include "weighting.fh"
#include "sbs.fh"
#include "print.fh"
#include "Molcas.fh"
#include "warnings.fh"
      Parameter(NRHS=1)
      Real*8 dSS(nInter,NRHS), rInt(nInter), dMass(nAtom), cMass(3),
     &       DFC(3*nAtom,NRHS), Coor(3,nAtom), Degen(3*nAtom)
      Character Lbl(nInter)*8, Name(nAtom)*(LENIN)
      Integer   iSym(3), iANr(nAtom),
     &          nStab(nAtom), jStab(0:7,nAtom), iCoSet(0:7,nAtom)
      Logical Smmtrc(3,nAtom), BSet, HSet, User_Def,
     &        Curvilinear, Numerical, DDV_Schlegel, Redundant,
     &        HWRS, Analytic_Hessian, PrQ, lOld, Invar, Error
      Save        BSet, HSet, lOld
*
      Data                 BSet/.False./, HSet/.False./,
     &     lOld/.False./
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
         Write (Lu,*) ' In Int2Car: Shifts'
         Write (Lu,*)
         Write (Lu,'(1X,A,2X,F10.4)') (Lbl(iInter),dss(iInter,1),
     &          iInter=1,nInter)
         Call RecPrt(' In Int2Car: qInt',' ',Work(ip_qInt),nInter,
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
         If (Abs(dss(i,1)) .gt. Abs(rMax)) Then
            rMax = dss(i,1)
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
      Do jter = 1, iterMx
*
*--------Compute the Cartesian shift, solve dq = B^T dx!
*
         M = 3*nAtom
         N = nInter
         Call Eq_Solver('T',M,N,NRHS,Work(ipBMx),Curvilinear,Degen,dSS,
     &                  DFC)
         Call Free_Work(ipBmx)
*
         If (iPrint.ge.99)
     &     Call PrList('Symmetry Distinct Nuclear Displacements',
     &                 Name,nAtom,DFC,3,nAtom)
*
*        Compute the RMS in Cartesian Coordinates.
*
         dx2=Zero
         denom=Zero
         Do ix = 1, 3*nAtom
            dx2 = dx2 + Degen(ix)*DFC(ix,1)**2
            denom=denom+Degen(ix)
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
         Call CofMss(Coor,dMass,nAtom,.False.,cMass,iSym)
         call dcopy_(3*nAtom,Coor,1,Cx(:,:,Iter+1),1)
         If (iPrint.ge.99)
     &      Call PrList('Symmetry Distinct Nuclear Coordinates / Bohr',
     &                   Name,nAtom,Coor,3,nAtom)
*
*--------Compute new values q and the Wilson B-matrix for the new
*        geometry with the current new set of Cartesian coordinates.
*
         nFix=0
         nWndw=1
         Call BMtrx(nLines,nBVct,ipBMx,nAtom,nInter,
     &              ip_qInt,Lbl,Coor,nDim,dMass,
     &              Name,Smmtrc,
     &              Degen,BSet,HSet,iter+1,ip_dqInt,
     &              Gx,mTtAtm,iANr,iOptH,User_Def,
     &              nStab,jStab,Curvilinear,Numerical,
     &              DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,
     &              PrQ,mxdc,iCoSet,lOld,rHidden,
     &              nFix,nQQ,iRef,Redundant,nqInt,MaxItr,nWndw)
*
*--------Check if the final structure is reached and get the
*        difference between the present structure and the final.
*
         rOld = rMax
         iMax_Old = iMax
         iMax = 1
         ip = ip_qInt + Iter*nInter -1
         rMax = Zero
         Do i = 1, nInter
            dSS(i,1) = rInt(i)-Work(ip+i)
            If (Abs(dSS(i,1)) .gt. Abs(rMax)) Then
               rMax = dSS(i,1)
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
            Write (Lu,'(1X,A,2X,F10.4)') (Lbl(iInter),dss(iInter,1),
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
      If (Error) Then
         Return
      Else
         Call WarningMessage(2,'Error in Int2Car')
         Write (Lu,*)
         Write (Lu,*) '***********************************************'
         Write (Lu,*) ' ERROR: No convergence in Int2Car !            '
         Write (Lu,*) ' Strong linear dependency among Coordinates.   '
         Write (Lu,*) ' Hint: Try to change the Internal Coordinates. '
         Write (Lu,*) '***********************************************'
         If (.NOT.User_Def) Call RecPrt('Int2Car: rInt  ','(10F15.10)',
     &                                  rInt,nInter,1)
         ip = ip_qInt + Iter*nInter
         If (.NOT.User_Def) Call RecPrt('Int2Car: qInt','(10F15.10)',
     &                                  Work(ip),nInter,1)
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
     &   Call Align(Cx(:,:,iter+1),Work(ipRef),nAtom)
*                                                                      *
************************************************************************
*                                                                      *
      Error=.False.
      Return
      End
