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
      SubRoutine Dissoc(xyz,nCntr,mCntr,rMss,nsAtom,Dist,B,lWrite,
     &                  Label,dB,ldB)
      Implicit Real*8 (A-H,O-Z)
************************************************************************
*                                                                      *
*     Object: To evaluate the B matrix elements of an internal         *
*             coordinate corresponding to a dissociation of two        *
*             parts of a molecule.                                     *
*                                                                      *
************************************************************************
#include "real.fh"
      Real*8 rMss(nCntr+mCntr), B(3,nCntr+mCntr), xyz(3,nCntr+mCntr),
     &       dB(3,nCntr+mCntr,3,nCntr+mCntr)
      Logical lWrite, ldB
      Real*8 R(3,2), RM(2)
      Character*8 Label
#include "angstr.fh"
*
      Call qEnter('Dissoc')
*
      call dcopy_(2,[Zero],0,RM,1)
      call dcopy_(6,[Zero],0,R,1)
*
#ifdef _DEBUGPRINT_
      Write (6,*) ' nCntr,mCntr=',nCntr,mCntr
      Call RecPrt(' Masses',' ',rMss,nCntr+mCntr,1)
      Call RecPrt(' xyz',' ',xyz,3,nCntr+mCntr)
#endif
      Do iCntr = 1, nCntr+mCntr
         i=1
         If (iCntr.gt.nCntr) i=2
*        Sum up mass of fragment
         RM(i) = RM(i) + rMss(iCntr)
*        Compute center of mass of the fragment
         Do ix = 1, 3
            R(ix,i) = R(ix,i) + rMss(iCntr) * xyz(ix,iCntr)
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('RM',' ',RM,1,2)
#endif
*
*     Evaluate center of mass of the two parts and the distance between
*
      Dist = Zero
      Do ix=1,3
         R(ix,1) = R(ix,1)/RM(1)
         R(ix,2) = R(ix,2)/RM(2)
         Dist = Dist + (R(ix,1)-R(ix,2))**2
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt(' Center of mass of fragments',' ',R,3,2)
#endif
*
      Dist = Sqrt(Dist)
*
      If (lWrite) Write (6,'(1X,A,A,2(F10.6,A))') Label,
     &   ' : Dissociation distance=',Dist,'/bohr',
     &       Dist*angstr,'/Angstrom'
*
*     Compute the B-matrix
*
      Do iCntr = 1, nCntr+mCntr
         If (iCntr.le.nCntr) Then
            Sign=1.0D0
            i=1
         Else
            Sign=-1.0D0
            i=2
         End If
         Do ix = 1, 3
            If (xyz(ix,iCntr).ne.Zero) Then
               Fact=Sign*rMss(iCntr)/RM(i)
            Else
               Fact=Zero
            End If
            B(ix,iCntr)= Fact* (R(ix,1)-R(ix,2)) / Dist
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('B',' ',B,3,nCntr+mCntr)
#endif
*
*     Compute the Cartesian derivative of the B-Matrix
*
      If (ldB) Then
         Call FZero(dB,(3*(nCntr+mCntr))**2)
         Do iCntr = 1, nCntr+mCntr
            If (iCntr.le.nCntr) Then
               Signi=1.0D0
               i=1
            Else
               Signi=-1.0D0
               i=2
            End If
            Facti=Signi*rMss(iCntr)/RM(i)
*
            Do jCntr = 1, nCntr+mCntr
               If (jCntr.le.nCntr) Then
                  Signj=1.0D0
                  j=1
               Else
                  Signj=-1.0D0
                  j=2
               End If
               Factj=Signj*rMss(jCntr)/RM(j)
*
               Do ix = 1, 3
                  If (xyz(ix,iCntr).ne.Zero) Then
                     dRdri=Facti
                  Else
                     dRdri=Zero
                  End If
                  Do jx = 1, 3
                     If (xyz(jx,jCntr).ne.Zero) Then
                        dRdrj=Factj
                     Else
                        dRdrj=Zero
                     End If
*
                     If (ix.eq.jx) Then
                        dB(ix,iCntr,jx,jCntr)=(dRdri*dRdrj
     &                                    -B(ix,iCntr)*B(jx,jCntr))/Dist
                     Else
                        dB(ix,iCntr,jx,jCntr)=(
     &                                    -B(ix,iCntr)*B(jx,jCntr))/Dist
                     End If
                  End Do
               End Do
*
            End Do
         End Do
#ifdef _DEBUGPRINT_
         Call RecPrt('dB',' ',dB,3*(nCntr+mCntr),3*(nCntr+mCntr))
#endif
      End If
      Call qExit('Dissoc')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nsAtom)
      End
