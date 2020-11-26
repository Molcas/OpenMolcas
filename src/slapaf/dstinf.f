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
      Subroutine DstInf(iStop,Just_Frequencies,Numerical)
      use Symmetry_Info, only: nIrrep, iOper
      use Slapaf_Info, only: Cx, Gx, Gx0, Coor, Weights
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "info_slapaf.fh"
#include "print.fh"
#include "SysDef.fh"
      Character*2 Element(MxAtom)
      Character*100 Get_SuperName, SuperName
      External Get_SuperName
#include "angstr.fh"
#include "weighting.fh"
      Integer, Allocatable:: Information(:)
      Real*8, Allocatable:: Cx_p(:,:), CC(:,:), RV(:,:), GxFix(:,:),
     &                      xyz(:,:)
*
      LOGICAL do_printcoords, do_fullprintcoords, Just_Frequencies,
     &        Found, Numerical
      Character*(LENIN) LblTMP(mxdc*nIrrep)
*                                                                      *
************************************************************************
*                                                                      *
*
      iRout=53
      iPrint=nPrint(iRout)
      do_printcoords=iPrint.ge.5
      do_fullprintcoords=(iPrint.gt.5.OR.iStop.gt.1)
      LOut=6
*                                                                      *
************************************************************************
*                                                                      *
*---  Write information of this iteration to the RLXITR file
      Call mma_allocate(Information,7,Label='Information')
      If (Stop) Then
         Information(1)=-99     ! Deactivate the record
         iOff_Iter=0
         Call Put_iScalar('iOff_Iter',iOff_Iter)
*
*        Restore the runfile data as if the computation was analytic
*        (note the gradient sign must be changed back)
*
         If (Just_Frequencies) Then
            Call Put_dScalar('Last Energy',Work(ipEner))
            Call mma_allocate(GxFix,3,nsAtom,Label='GxFix')
            call dcopy_(3*nsAtom,Gx,1,GxFix,1)
            GxFix(:,:) = - GxFix(:,:)
            Call Put_Grad(GxFix,3*nsAtom)
            Call mma_deallocate(GxFix)
            Call Put_dArray('Unique Coordinates',Cx,3*nsAtom)
            Call Put_Coord_New(Cx,nsAtom)
         End If
      Else
         Call qpg_iArray('Slapaf Info 1',Found,nSlap)
         If (Found) Then
            Call Get_iArray('Slapaf Info 1',Information,7)
            If (Information(1).ne.-99) Information(1)=MaxItr
         Else
            Information(1)=MaxItr
         End If
      End If
*
      SuperName=Get_Supername()
      If (SuperName.ne.'numerical_gradient') Then
         Information(2)=Iter
         Information(3)=mTROld ! # symm. transl /rot.
         If (lOld_Implicit) Then
            Information(4)=1
         Else
            Information(4)=0
         End If
         Information(5)=ipEner-ipRlx
         Information(6)=ipCx-ipRlx
         Information(7)=ipGx-ipRlx
         Call Put_iArray('Slapaf Info 1',Information,7)
         Call mma_deallocate(Information)

         Call Dcopy_(3*nsAtom*(MaxItr+1),Cx,1,Work(ipCx),1)
         Call Dcopy_(3*nsAtom*(MaxItr+1),Gx,1,Work(ipGx),1)
         Call Dcopy_(3*nsAtom*(MaxItr+1),Gx0,1,Work(ipGx0),1)

         Call Put_dArray('Slapaf Info 2',Work(ipRlx),Lngth)
         Call Put_cArray('Slapaf Info 3',Stat(0),(MaxItr+1)*128)
         Call Put_dArray('qInt',Work(ipqInt),nqInt)
         Call Put_dArray('dqInt',Work(ipdqInt),nqInt)
      End If
      If (Just_Frequencies) Return
*                                                                      *
************************************************************************
*                                                                      *
*---- Geometry information
*
      If (Stop.or.do_printcoords) Then
         Write (LOut,*)
         Call CollapseOutput(1,'Geometry section')
         Write (LOut,*)
         Write (LOut,'(80A)') ('*',i=1,80)
         If (Stop) Then
            Write (LOut,*)
     &        ' Geometrical information of the final structure'
            r_Iter=DBLE(Iter)
            Call Add_Info('GEO_ITER',[r_Iter],1,8)
         Else IF (do_printcoords) THEN
            Write (LOut,*)
     &         ' Geometrical information of the new structure'
         End If
         Write (LOut,'(80A)') ('*',i=1,80)
         Write (LOut,*)
      End If
*
      Call Get_iScalar('Pseudo atoms',nsAtom_p)
      If (nsAtom_p.gt.0) Then
         Call mma_allocate(Cx_p,3,nsAtom_p,Label='Cx_p')
         Call Get_dArray('Pseudo Coordinates',Cx_p,3*nsAtom_p)
      End If
*
      Call mma_Allocate(CC,3,nIrrep*(nsAtom+nsAtom_p),Label='CC')
      nTemp = 0
      Do isAtom = 1, nsAtom + nsAtom_p
         If (isAtom.le.nsAtom) Then
            x1 = Coor(1,isAtom)
            y1 = Coor(2,isAtom)
            z1 = Coor(3,isAtom)
         Else
            jsAtom = isAtom - nsAtom
            x1 = Cx_p(1,jsAtom)
            y1 = Cx_p(2,jsAtom)
            z1 = Cx_p(3,jsAtom)
         End If
         Do 6001 iIrrep = 0, nIrrep-1
            x2 = x1
            If (iAnd(1,iOper(iIrrep)).ne.0) x2 = - x2
            y2 = y1
            If (iAnd(2,iOper(iIrrep)).ne.0) y2 = - y2
            z2 = z1
            If (iAnd(4,iOper(iIrrep)).ne.0) z2 = - z2
*
*           Check if it is already in the list
*
            Do iTemp = 1, nTemp
               r = (x2-CC(1,iTemp))**2 +
     &             (y2-CC(2,iTemp))**2 +
     &             (z2-CC(3,iTemp))**2
               If (r.eq.Zero) Go To 6001
            End Do
            nTemp = nTemp + 1
            If (nTemp.gt.nIrrep*(nsAtom+nsAtom_p)) Then
               Call WarningMessage(2,'Error in DstInf')
               Write (6,*) 'nTemp.gt.nIrrep*nsAtom'
               Call Abend()
            End If
            CC(1,nTemp) = x2
            CC(2,nTemp) = y2
            CC(3,nTemp) = z2
            If (isAtom.le.nsAtom) Then
               LblTMP(nTemp)=AtomLbl(isAtom)
            Else
               LblTMP(nTemp)='PC'
            End If
 6001    Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*-----Write out the new cartesian symmetry coordinates.
*
      If (Stop) Then
         Write (LOut,*) ' NOTE: on convergence the final predicted '
     &               //'structure will be printed here.'
         Write (LOut,*) ' This is not identical to the structure'
     &               //' printed in the head of the output.'
         Call OutCoor(
     &    '* Nuclear coordinates of the final structure / Bohr     *',
     &    AtomLbl,nsAtom,Coor,3,nsAtom,.False.)
         Call OutCoor(
     &    '* Nuclear coordinates of the final structure / Angstrom *',
     &    AtomLbl,nsAtom,Coor,3,nsAtom,.True.)
      Else If (Do_PrintCoords) Then
         Call OutCoor(
     &    '* Nuclear coordinates for the next iteration / Bohr     *',
     &    AtomLbl,nsAtom,Coor,3,nsAtom,.False.)
         Call OutCoor(
     &    '* Nuclear coordinates for the next iteration / Angstrom *',
     &    AtomLbl,nsAtom,Coor,3,nsAtom,.True.)
      End If
*
      If (nsAtom_p.gt.0) Then
         iOff = nTemp - nsAtom_p + 1
         Call OutCoor(
     &'* Pseudo charge coordinates for the next iteration / Bohr     *',
     &    LblTMP(iOff),nsAtom_p,Cx_p,3,nsAtom_p,.False.)
         Call OutCoor(
     &'* Pseudo Charge coordinates for the next iteration / Angstrom *',
     &    LblTMP(iOff),nsAtom_p,Cx_p,3,nsAtom_p,.True.)
         Call mma_deallocate(Cx_p)
      End If
*
      IF (do_printcoords) THEN
         Call Get_iScalar('N ZMAT',N_ZMAT)
         If (N_ZMAT.GT.0) Call OutZMAT(nsAtom,Coor,N_ZMAT)
*
         IF (do_fullprintcoords) THEN
           If (nTemp.ge.2)
     &     Call Dstncs(LblTMP,CC,nTemp,angstr,Max_Center,5)
*
           If (nTemp.ge.3)
     &     Call Angles(LblTMP,CC,nTemp,Rtrnc,Max_Center)
*
           If (nTemp.ge.4)
     &     Call Dihedr(LblTMP,CC,nTemp,Rtrnc,Max_Center)
         END IF
      END IF
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(CC)
*                                                                      *
************************************************************************
*                                                                      *
*---  Write the new cartesian symmetry coordinates on GEONEW
*
      Call Put_Coord_New(Cx(1,1,iter+1),nsAtom)
*                                                                      *
************************************************************************
*                                                                      *
*     If two runfiles are associated with the calculation update both
*     files.
*
      Call f_Inquire('RUNFILE2',Found)
      If (Found) Then
         Call NameRun('RUNFILE2')
         Call Put_Coord_New(Cx(1,1,iter+1),nsAtom)
         Call NameRun('RUNFILE')
      End If
*
*     Update the .Opt.xyz file
*
      If (.Not.Numerical) Then
         Call Get_nAtoms_All(nCoord)
         Call mma_allocate(xyz,3,nCoord,Label='xyz')
         Call Get_Coord_New_All(xyz,nCoord)
         Call Get_Name_All(Element)
*
         Lu_xyz=IsFreeUnit(11)
         Call MOLCAS_Open(Lu_xyz,'XYZ')
         Write (Lu_xyz,'(I4)') nCoord
         Write(Lu_xyz,*) Work(ipEner+Iter)
*        Write (Lu_xyz,'(A)') 'Coordinates generated by Slapaf'
         Do i = 1, nCoord
            Write (Lu_xyz,'(A2,3F15.8)') Element(i),
     &            (Angstr*xyz(j,i),j=1,3)
         End Do
         Close (Lu_xyz)
         Call mma_deallocate(xyz)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     If a transition state optimization put the "reaction" vector
*     on the RUNFILE(S)
*
       If (iAnd(iOptC,128).ne.128 .and. Stop ) Then
*
           Call mma_allocate(RV,3,nsAtom,Label='RV')
           Call dcopy_(3*nsAtom,Work(ipMF),1,RV,1)
           Do i=1,nsAtom
             xWeight=Weights(i)
             RV(:,i) = RV(:,i)/xWeight
           End Do
           Call OutCoor('* The Cartesian Reaction vector'//
     &                  '                         *',
     &                  AtomLbl,nsAtom,RV,3,nsAtom,.True.)

           Call f_Inquire('RUNREAC',Found)
           If (Found) Then
              Call NameRun('RUNREAC')
              Call Put_dArray('Reaction Vector',RV,3*nsAtom)
           End If
           Call f_Inquire('RUNPROD',Found)
           If (Found) Then
              Call NameRun('RUNPROD')
              Call Put_dArray('Reaction Vector',RV,3*nsAtom)
           End If
           Call NameRun('RUNFILE')
           Call Put_dArray('Reaction Vector',RV,3*nsAtom)
           Call mma_deallocate(RV)
           iDo_dDipM=0
           Call GF_on_the_fly(iDo_dDipM)
*
       End If
*                                                                      *
************************************************************************
*                                                                      *
      If (Stop.or.do_printcoords)
     &   Call CollapseOutput(0,'Geometry section')
      Return
      End
