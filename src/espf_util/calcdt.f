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
      Subroutine CalcDT(nMult,nGrdPt,natom,nAtQM,ipIsMM,iGrdTyp,Coord,
     &                  Grid,DGrid,T,TT,DT,DTT,DTTTT,DTTT)
      Implicit Real*8 (A-H,O-Z)
*
*     The dependancy of T with respect to grid points is taken
*     into account if the GEPOL grid is used (iGrdTyp == 2)
*
*     1) dT    = dT/dx + Sum_k dT/dq_k * dq_k/dx
*     2) dTT   = (dTt.T+Tt.dT)
*     3) dTTTT = (dTt.T+Tt.dT)((TtT)**-1)
*     4) dTT   = d((TtT)**-1) = -((TtT)**-1)(dTt.T+Tt.dT)((TtT)**-1)
*     5) dTTT = d((TtT)**-1)*Tt + ((TtT)**-1)*dTt
*
*
#include "espf.fh"
*
      Dimension Coord(3,natom),Grid(3,nGrdPt),DGrid(nGrdPt,nAtQM,3,3),
     &          T(nMult,nGrdPt),TT(nMult,nMult),
     &          DT(nMult,nGrdPt,3,nAtQM),DTT(nMult,nMult,3,nAtQM),
     &          DTTTT(nMult,nMult,3,nAtQM),DTTT(nMult,nGrdPt,3,nAtQM)
*
*     Very local array
*
      Dimension DG(3,3)
*
      Call QEnter('calcdt')
      iPL = iPL_espf()
*
      nOrd = nMult/nAtQM
*
*     Step 1
*
      iQM = 0
      Do iPnt = 1, nGrdPt
         iQM = 0
         Do iAt = 1, natom
            If (iWork(ipIsMM+iAt-1).eq.1) Goto 1
            iQM = iQM + 1
            X = Grid(1,iPnt) - Coord(1,iAt)
            Y = Grid(2,iPnt) - Coord(2,iAt)
            Z = Grid(3,iPnt) - Coord(3,iAt)
            R = sqrt(X*X+Y*Y+Z*Z)
            R2 = R * R
            R3 = R2 * R
            R5 = R2 * R3
            Do jQM = 1, nAtQM
               Do I = 1, 3
                  Do J = 1, 3
                     DG(I,J) = Zero
*
* This is commented since the derivative of the GEPOL grid points look ugly
*
c                     If (iGrdTyp.eq.2) DG(I,J) = DGrid(iPnt,jQM,I,J)
                  End Do
               End Do
               dIJ = Zero
               If (iQM .eq. jQM) dIJ = One
               If (nOrd .eq. 1) Then
                  DT(iQM   ,iPnt,1,jQM) = X/R3*(dIJ-DG(1,1))
     &                                  + Y/R3*(   -DG(2,1))
     &                                  + Z/R3*(   -DG(3,1))
                  DT(iQM   ,iPnt,2,jQM) = X/R3*(   -DG(1,2))
     &                                  + Y/R3*(dIJ-DG(2,2))
     &                                  + Z/R3*(   -DG(3,2))
                  DT(iQM   ,iPnt,3,jQM) = X/R3*(   -DG(1,3))
     &                                  + Y/R3*(   -DG(2,3))
     &                                  + Z/R3*(dIJ-DG(3,3))
               Else
                 DT(4*iQM-3,iPnt,1,jQM)= X/R3*(dIJ-DG(1,1))
     &                                 + Y/R3*(   -DG(2,1))
     &                                 + Z/R3*(   -DG(3,1))
                 DT(4*iQM-3,iPnt,2,jQM)= X/R3*(   -DG(1,2))
     &                                 + Y/R3*(dIJ-DG(2,2))
     &                                 + Z/R3*(   -DG(3,2))
                 DT(4*iQM-3,iPnt,3,jQM)= X/R3*(   -DG(1,3))
     &                                 + Y/R3*(   -DG(2,3))
     &                                 + Z/R3*(dIJ-DG(3,3))
                 DT(4*iQM-2,iPnt,1,jQM)= (three*X*X-R2)/R5*(dIJ-DG(1,1))
     &                                 + (three*X*Y   )/R5*(   -DG(2,1))
     &                                 + (three*X*Z   )/R5*(   -DG(3,1))
                 DT(4*iQM-2,iPnt,2,jQM)= (three*X*X-R2)/R5*(   -DG(1,2))
     &                                 + (three*X*Y   )/R5*(dIJ-DG(2,2))
     &                                 + (three*X*Z   )/R5*(   -DG(3,2))
                 DT(4*iQM-2,iPnt,3,jQM)= (three*X*X-R2)/R5*(   -DG(1,3))
     &                                 + (three*X*Y   )/R5*(   -DG(2,3))
     &                                 + (three*X*Z   )/R5*(dIJ-DG(3,3))
                 DT(4*iQM-1,iPnt,1,jQM)= (three*Y*X   )/R5*(dIJ-DG(1,1))
     &                                 + (three*Y*Y-R2)/R5*(   -DG(2,1))
     &                                 + (three*Y*Z   )/R5*(   -DG(3,1))
                 DT(4*iQM-1,iPnt,2,jQM)= (three*Y*X   )/R5*(   -DG(1,2))
     &                                 + (three*Y*Y-R2)/R5*(dIJ-DG(2,2))
     &                                 + (three*Y*Z   )/R5*(   -DG(3,2))
                 DT(4*iQM-1,iPnt,3,jQM)= (three*Y*X   )/R5*(   -DG(1,3))
     &                                 + (three*Y*Y-R2)/R5*(   -DG(2,3))
     &                                 + (three*Y*Z   )/R5*(dIJ-DG(3,3))
                 DT(4*iQM  ,iPnt,1,jQM)= (three*Z*X   )/R5*(dIJ-DG(1,1))
     &                                 + (three*Z*Y   )/R5*(   -DG(2,1))
     &                                 + (three*Z*Z-R2)/R5*(   -DG(3,1))
                 DT(4*iQM  ,iPnt,2,jQM)= (three*Z*X   )/R5*(   -DG(1,2))
     &                                 + (three*Z*Y   )/R5*(dIJ-DG(2,2))
     &                                 + (three*Z*Z-R2)/R5*(   -DG(3,2))
                 DT(4*iQM  ,iPnt,3,jQM)= (three*Z*X   )/R5*(   -DG(1,3))
     &                                 + (three*Z*Y   )/R5*(   -DG(2,3))
     &                                 + (three*Z*Z-R2)/R5*(dIJ-DG(3,3))
               EndIf
            End Do
1           Continue
         End Do
      End Do
      If (iQM.ne.nAtQM) Then
         Write(6,'(A)') ' Error in espf/initdb: iQM != nAtQM !!!'
         Call Abend()
      End If
      If (iPL.ge.4) Then
         Do iMlt = 1, nMult
            Do jPnt = 1, nGrdPt
               Write(6,'(A,i4,i4)') ' DT for iMlt and jPnt: ',iMlt,jPnt
               Do iQM = 1, nAtQM
                  Write(6,'(i4,3F20.7)') iQM,(DT(iMlt,jPnt,I,iQM),I=1,3)
               End Do
            End Do
         End Do
      End If
*
*     Step 2
*
      Do iMlt = 1, nMult
         Do jMlt = 1, nMult
            Do iQM = 1, nAtQM
               DTT(iMlt,jMlt,1,iQM)= Zero
               DTT(iMlt,jMlt,2,iQM)= Zero
               DTT(iMlt,jMlt,3,iQM)= Zero
               Do iPnt = 1, nGrdPt
                  DTT(iMlt,jMlt,1,iQM)= DTT(iMlt,jMlt,1,iQM)
     &                                + T(iMlt,iPnt)*DT(jMlt,iPnt,1,iQM)
     &                                + DT(iMlt,iPnt,1,iQM)*T(jMlt,iPnt)
                  DTT(iMlt,jMlt,2,iQM)= DTT(iMlt,jMlt,2,iQM)
     &                                + T(iMlt,iPnt)*DT(jMlt,iPnt,2,iQM)
     &                                + DT(iMlt,iPnt,2,iQM)*T(jMlt,iPnt)
                  DTT(iMlt,jMlt,3,iQM)= DTT(iMlt,jMlt,3,iQM)
     &                                + T(iMlt,iPnt)*DT(jMlt,iPnt,3,iQM)
     &                                + DT(iMlt,iPnt,3,iQM)*T(jMlt,iPnt)
               End Do
            End Do
            If (iPL.ge.4) Then
               Write(6,'(A,i4,i4)') 'DTT1 for iMlt and jMlt: ',iMlt,jMlt
               Do iQM = 1, nAtQM
                 Write(6,'(i4,3F20.7)') iQM,(DTT(iMlt,jMlt,I,iQM),I=1,3)
               End Do
            End If
         End Do
      End Do
*
*     Step 3
*
      Do iMlt = 1, nMult
         Do jMlt = 1, nMult
            Do iQM = 1, nAtQM
               DTTTT(iMlt,jMlt,1,iQM) = Zero
               DTTTT(iMlt,jMlt,2,iQM) = Zero
               DTTTT(iMlt,jMlt,3,iQM) = Zero
               Do kMlt = 1, nMult
                  DTTTT(iMlt,jMlt,1,iQM)= DTTTT(iMlt,jMlt,1,iQM)
     &                              + DTT(iMlt,kMlt,1,iQM)*TT(kMlt,jMlt)
                  DTTTT(iMlt,jMlt,2,iQM)= DTTTT(iMlt,jMlt,2,iQM)
     &                              + DTT(iMlt,kMlt,2,iQM)*TT(kMlt,jMlt)
                  DTTTT(iMlt,jMlt,3,iQM)= DTTTT(iMlt,jMlt,3,iQM)
     &                              + DTT(iMlt,kMlt,3,iQM)*TT(kMlt,jMlt)
               End Do
            End Do
            If (iPL.ge.4) Then
              Write(6,'(A,i4,i4)') 'DTTTT for iMlt and jMlt: ',iMlt,jMlt
              Do iQM = 1, nAtQM
               Write(6,'(i4,3F20.7)') iQM,(DTTTT(iMlt,jMlt,I,iQM),I=1,3)
              End Do
           End If
         End Do
      End Do
*
*     Step 4
*
      Do iMlt = 1, nMult
         Do jMlt = 1, nMult
            Do iQM = 1, nAtQM
               DTT(iMlt,jMlt,1,iQM) = Zero
               DTT(iMlt,jMlt,2,iQM) = Zero
               DTT(iMlt,jMlt,3,iQM) = Zero
               Do kMlt = 1, nMult
                  DTT(iMlt,jMlt,1,iQM)= DTT(iMlt,jMlt,1,iQM)
     &                            - TT(iMlt,kMlt)*DTTTT(kMlt,jMlt,1,iQM)
                  DTT(iMlt,jMlt,2,iQM)= DTT(iMlt,jMlt,2,iQM)
     &                            - TT(iMlt,kMlt)*DTTTT(kMlt,jMlt,2,iQM)
                  DTT(iMlt,jMlt,3,iQM)= DTT(iMlt,jMlt,3,iQM)
     &                            - TT(iMlt,kMlt)*DTTTT(kMlt,jMlt,3,iQM)
               End Do
            End Do
            If (iPL.ge.4) Then
               Write(6,'(A,i4,i4)') 'DTT for iMlt and jMlt: ',iMlt,jMlt
               Do iQM = 1, nAtQM
                 Write(6,'(i4,3F20.7)') iQM,(DTT(iMlt,jMlt,I,iQM),I=1,3)
               End Do
            End If
         End Do
      End Do
*
*     Step 5
*
      Do iMlt = 1, nMult
         Do iPnt = 1, nGrdPt
            Do iQM = 1, nAtQM
               DTTT(iMlt,iPnt,1,iQM) = Zero
               DTTT(iMlt,iPnt,2,iQM) = Zero
               DTTT(iMlt,iPnt,3,iQM) = Zero
               Do jMlt = 1, nMult
                  DTTT(iMlt,iPnt,1,iQM)= DTTT(iMlt,iPnt,1,iQM)
     &                               + TT(iMlt,jMlt)*DT(jMlt,iPnt,1,iQM)
     &                               + DTT(iMlt,jMlt,1,iQM)*T(jMlt,iPnt)
                  DTTT(iMlt,iPnt,2,iQM)= DTTT(iMlt,iPnt,2,iQM)
     &                               + TT(iMlt,jMlt)*DT(jMlt,iPnt,2,iQM)
     &                               + DTT(iMlt,jMlt,2,iQM)*T(jMlt,iPnt)
                  DTTT(iMlt,iPnt,3,iQM)= DTTT(iMlt,iPnt,3,iQM)
     &                               + TT(iMlt,jMlt)*DT(jMlt,iPnt,3,iQM)
     &                               + DTT(iMlt,jMlt,3,iQM)*T(jMlt,iPnt)
               End Do
            End Do
            If (iPL.ge.4) Then
               Write(6,'(A,i4,i4)') 'DTTT for iMlt and iPnt: ',iMlt,iPnt
               Do iQM = 1, nAtQM
                Write(6,'(i4,3F20.7)') iQM,(DTTT(iMlt,iPnt,I,iQM),I=1,3)
               End Do
            End If
         End Do
      End Do
*
      Call QExit('calcdt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(iGrdTyp)
        Call Unused_real_array(DGrid)
      End If
      End
