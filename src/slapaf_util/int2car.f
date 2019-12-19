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
     &                  nDim,Lbl,Name,iOper,nSym,iSym,
     &                  Smmtrc,Degen,iter,
     &                  ip_dqInt,Gx,Cx,mTtAtm,iANr,iOptH,
     &                  User_Def,nStab,jStab,Curvilinear,
     &                  Numerical,DDV_Schlegel,HWRS,
     &                  Analytic_Hessian,iOptC,PrQ,mxdc,
     &                  iCoSet,rHidden,Error,ipRef,Redundant,nqInt,
     &                  MaxItr)
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
#include "print.fh"
#include "Molcas.fh"
#include "warnings.fh"
      Parameter(NRHS=1)
      Real*8 dSS(nInter,NRHS), rInt(nInter), dMass(nAtom),
     &       DFC(3*nAtom,NRHS),
     &       Coor(3,nAtom), Degen(3*nAtom),
     &       Gx(3*nAtom,iter), Cx(3*nAtom,iter+1),
     &       cMass(3)
      Character Lbl(nInter)*8, Name(nAtom)*(LENIN)
      Integer   iOper(0:7), iSym(3), iANr(nAtom),
     &          nStab(nAtom), jStab(0:7,nAtom), iCoSet(0:7,nAtom)
      Logical Smmtrc(3,nAtom), BSet, HSet, User_Def,
     &        Curvilinear, Numerical, DDV_Schlegel, Redundant,
     &        HWRS, Analytic_Hessian, PrQ, lOld
      Save        BSet, HSet, lOld
*
C     Data Error/1.0D-06/, BSet/.False./, HSet/.False./,
      Data                 BSet/.False./, HSet/.False./,
     &     lOld/.False./
*                                                                      *
************************************************************************
*                                                                      *
      Lu=6
      iRout = 33
      iPrint=nPrint(iRout)
      Call QEnter('Int2Car')
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
            If (Cx((iAtom-1)*3+1,Iter).eq.Zero .and.
     &          Abs(Coor(1,iAtom)).lt.1.0D-13) Coor(1,iAtom)=Zero
            If (Cx((iAtom-1)*3+2,Iter).eq.Zero .and.
     &          Abs(Coor(2,iAtom)).lt.1.0D-13) Coor(2,iAtom)=Zero
            If (Cx((iAtom-1)*3+3,Iter).eq.Zero .and.
     &          Abs(Coor(3,iAtom)).lt.1.0D-13) Coor(3,iAtom)=Zero
         End Do
*
         Call CofMss(Coor,dMass,iOper,nSym,nAtom,.False.,cMass,iSym)
         call dcopy_(3*nAtom,Coor,1,Cx(1,Iter+1),1)
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
     &              Name,nSym,iOper,Smmtrc,
     &              Degen,BSet,HSet,iter+1,ip_dqInt,
     &              ipShift,Gx,Cx,mTtAtm,iANr,iOptH,User_Def,
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
      Call WarningMessage(2,'Error in Int2Car')
      Write (Lu,*)
      Write (Lu,*) '***********************************************'
      Write (Lu,*) ' ERROR: No convergence in Int2Car !            '
      Write (Lu,*) ' Strong linear dependency among Coordinates.   '
      Write (Lu,*) ' Hint: Try to change the Internal Coordinates. '
      Write (Lu,*) '***********************************************'
      If (.NOT.User_Def) Call RecPrt('Int2Car: rInt  ','(10F15.10)',
     &                               rInt,nInter,1)
      ip = ip_qInt + Iter*nInter
      If (.NOT.User_Def) Call RecPrt('Int2Car: qInt','(10F15.10)',
     &                               Work(ip),nInter,1)
      Write (Lu,*)
      Call Quit(_RC_NOT_CONVERGED_)
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
*     In case the internal coordinates where nonvalid we now take the
*     final internal coordinates at the converged Cartesian structure
*     as the internal coordinates.
*
      ip = ip_qInt + Iter*nInter
      call dcopy_(nInter,rInt,1,Work(ip),1)
*                                                                      *
************************************************************************
*                                                                      *
*     Finally, just to be safe align the new Cartesian structure with
*     the reference structure (see init2.f)
*
      If (WeightedConstraints)
     &   Call Align(Cx(1,iter+1),Work(ipRef),nAtom)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('Int2Car')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(Error)
      End
