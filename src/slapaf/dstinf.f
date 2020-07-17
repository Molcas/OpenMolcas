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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "info_slapaf.fh"
#include "print.fh"
#include "SysDef.fh"
      Character*2 Element(MxAtom)
      Character*100 Get_SuperName, SuperName
      External Get_SuperName
#include "angstr.fh"
#include "weighting.fh"
*
      LOGICAL do_printcoords, do_fullprintcoords, Just_Frequencies,
     &        Found, Numerical
      Character*(LENIN) LblTMP(mxdc*nSym)
*                                                                      *
************************************************************************
*                                                                      *
*     Call QEnter('DstInf')
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
      Call GetMem(' iter','Allo','Inte',ipItr,7)
      If (Stop) Then
         iWork(ipItr)=-99     ! Deactivate the record
         iOff_Iter=0
         Call Put_iScalar('iOff_Iter',iOff_Iter)
*
*        Restore the runfile data as if the computation was analytic
*        (note the gradient sign must be changed back)
*
         If (Just_Frequencies) Then
            Call Put_dScalar('Last Energy',Work(ipEner))
            Call Allocate_Work(ipGxFix,3*nsAtom)
            call dcopy_(3*nsAtom,Work(ipGx),1,Work(ipGxFix),1)
            Call DScal_(3*nsAtom,-One,Work(ipGxFix),1)
            Call Put_Grad(Work(ipGxFix),3*nsAtom)
            Call Free_Work(ipGxFix)
            Call Put_dArray('Unique Coordinates',Work(ipCx),3*nsAtom)
            Call Put_Coord_New(Work(ipCx),nsAtom)
         End If
      Else
         Call qpg_iArray('Slapaf Info 1',Found,nSlap)
         If (Found) Then
            Call Get_iArray('Slapaf Info 1',iWork(ipItr),7)
            If (iWork(ipItr).ne.-99) iWork(ipItr)=MaxItr
         Else
            iWork(ipItr)=MaxItr
         End If
      End If
*
      SuperName=Get_Supername()
      If (SuperName.ne.'numerical_gradient') Then
         iWork(ipItr+1)=Iter
         iWork(ipItr+2)=mTROld ! # symm. transl /rot.
         If (lOld_Implicit) Then
            iWork(ipItr+3)=1
         Else
            iWork(ipItr+3)=0
         End If
         iWork(ipItr+4)=ipEner-ipRlx
         iWork(ipItr+5)=ipCx-ipRlx
         iWork(ipItr+6)=ipGx-ipRlx
         Call Put_iArray('Slapaf Info 1',iWork(ipItr),7)
         Call GetMem(' iter','Free','Inte',ipItr,7)
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
         Call GetMem('Coor_p','Allo','Real',ipCx_p,3*nsAtom_p)
         Call Get_dArray('Pseudo Coordinates',Work(ipCx_p),3*nsAtom_p)
      Else
         ipCx_p=ip_Dummy
      End If
*
      Call GetMem('Carcor','Allo','Real',ipCC,3*nSym*(nsAtom+nsAtom_p))
      nTemp = 0
      Do isAtom = 1, nsAtom + nsAtom_p
         If (isAtom.le.nsAtom) Then
            x1 = Work(ipCoor-1+(isAtom-1)*3+1)
            y1 = Work(ipCoor-1+(isAtom-1)*3+2)
            z1 = Work(ipCoor-1+(isAtom-1)*3+3)
         Else
            x1 = Work(ipCx_p-1+(isAtom-1-nsAtom)*3+1)
            y1 = Work(ipCx_p-1+(isAtom-1-nsAtom)*3+2)
            z1 = Work(ipCx_p-1+(isAtom-1-nsAtom)*3+3)
         End If
         Do 6001 iIrrep = 0, nSym-1
            x2 = x1
            If (iAnd(1,iOper(iIrrep)).ne.0) x2 = - x2
            y2 = y1
            If (iAnd(2,iOper(iIrrep)).ne.0) y2 = - y2
            z2 = z1
            If (iAnd(4,iOper(iIrrep)).ne.0) z2 = - z2
            Do iTemp = 1, nTemp
               r = (x2-Work(ipCC-1+(iTemp-1)*3+1))**2 +
     &             (y2-Work(ipCC-1+(iTemp-1)*3+2))**2 +
     &             (z2-Work(ipCC-1+(iTemp-1)*3+3))**2
               If (r.eq.Zero) Go To 6001
            End Do
            nTemp = nTemp + 1
            If (nTemp.gt.nSym*(nsAtom+nsAtom_p)) Then
               Call WarningMessage(2,'Error in DstInf')
               Write (6,*) 'nTemp.gt.nSym*nsAtom'
               Call Abend()
            End If
            Work(ipCC-1+(nTemp-1)*3+1) = x2
            Work(ipCC-1+(nTemp-1)*3+2) = y2
            Work(ipCC-1+(nTemp-1)*3+3) = z2
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
     &    AtomLbl,nsAtom,Work(ipCoor),3,nsAtom,.False.)
         Call OutCoor(
     &    '* Nuclear coordinates of the final structure / Angstrom *',
     &    AtomLbl,nsAtom,Work(ipCoor),3,nsAtom,.True.)
      Else If (Do_PrintCoords) Then
         Call OutCoor(
     &    '* Nuclear coordinates for the next iteration / Bohr     *',
     &    AtomLbl,nsAtom,Work(ipCoor),3,nsAtom,.False.)
         Call OutCoor(
     &    '* Nuclear coordinates for the next iteration / Angstrom *',
     &    AtomLbl,nsAtom,Work(ipCoor),3,nsAtom,.True.)
      End If
*
      If (nsAtom_p.gt.0) Then
         iOff = nTemp - nsAtom_p + 1
         Call OutCoor(
     &'* Pseudo charge coordinates for the next iteration / Bohr     *',
     &    LblTMP(iOff),nsAtom_p,Work(ipCx_p),3,nsAtom_p,.False.)
         Call OutCoor(
     &'* Pseudo Charge coordinates for the next iteration / Angstrom *',
     &    LblTMP(iOff),nsAtom_p,Work(ipCx_p),3,nsAtom_p,.True.)
         Call Free_Work(ipCx_p)
      End If
*
      IF (do_printcoords) THEN
         Call Get_iScalar('N ZMAT',N_ZMAT)
         If (N_ZMAT.GT.0) Call OutZMAT(nsAtom,Work(ipCoor),N_ZMAT)
*
         IF (do_fullprintcoords) THEN
           If (nTemp.ge.2)
     &     Call Dstncs(LblTMP,Work(ipCC),nTemp,
     &                 angstr,Max_Center,5)
*
           If (nTemp.ge.3)
     &     Call Angles(LblTMP,Work(ipCC),nTemp,Rtrnc,Max_Center)
*
           If (nTemp.ge.4)
     &     Call Dihedr(LblTMP,Work(ipCC),nTemp,Rtrnc,Max_Center)
         END IF
      END IF
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('Carcor','Free','Real',ipCC,3*nSym*nsAtom)
*                                                                      *
************************************************************************
*                                                                      *
*---  Write the new cartesian symmetry coordinates on GEONEW
*
      jpCoor = ipCx + iter*3*nsAtom
      Call Put_Coord_New(Work(jpCoor),nsAtom)
      call recprt('CoordNew','',Work(jpCoor),nsAtom,3)
      call recprt('ipCoor','',Work(ipCoor),nsAtom,3)
*                                                                      *
************************************************************************
*                                                                      *
*     If two runfiles are associated with the calculation update both
*     files.
*
      Call f_Inquire('RUNFILE2',Found)
      If (Found) Then
         Call NameRun('RUNFILE2')
         Call Put_Coord_New(Work(jpCoor),nsAtom)
         Call NameRun('RUNFILE')
      End If
*
*     Update the .Opt.xyz file
*
      If (.Not.Numerical) Then
         Call Get_nAtoms_All(nCoord)
         Call Allocate_Work(ipxyz,3*nCoord)
         Call Get_Coord_New_All(Work(ipxyz),nCoord)
         Call Get_Name_All(Element)
*
         Lu_xyz=IsFreeUnit(11)
         Call MOLCAS_Open(Lu_xyz,'XYZ')
         Write (Lu_xyz,'(I4)') nCoord
         Write(Lu_xyz,*) Work(ipEner+Iter)
*        Write (Lu_xyz,'(A)') 'Coordinates generated by Slapaf'
         Do i = 1, nCoord
            Write (Lu_xyz,'(A2,3F15.8)') Element(i),
     &            (Angstr*Work(ipxyz+(i-1)*3+j),j=0,2)
         End Do
         Close (Lu_xyz)
      Call Free_Work(ipxyz)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     If a transition state optimization put the "reaction" vector
*     on the RUNFILE(S)
*
       If (iAnd(iOptC,128).ne.128 .and. Stop ) Then
*
           Call GetMem('ReacV','Allo','Real',ipRV,3*nsAtom)
           Call dcopy_(3*nsAtom,Work(ipMF),1,Work(ipRV),1)
           Do i=0,nsAtom-1
             xWeight=Work(ipWeights+i)
             Call DScal_(3,One/xWeight,Work(ipRV+3*i),1)
           End Do
           Call OutCoor('* The Cartesian Reaction vector'//
     &                  '                         *',
     &                  AtomLbl,nsAtom,Work(ipRV),3,nsAtom,.True.)

           Call f_Inquire('RUNREAC',Found)
           If (Found) Then
              Call NameRun('RUNREAC')
              Call Put_dArray('Reaction Vector',Work(ipRV),3*nsAtom)
           End If
           Call f_Inquire('RUNPROD',Found)
           If (Found) Then
              Call NameRun('RUNPROD')
              Call Put_dArray('Reaction Vector',Work(ipRV),3*nsAtom)
           End If
           Call NameRun('RUNFILE')
           Call Put_dArray('Reaction Vector',Work(ipRV),3*nsAtom)
           Call GetMem('ReacV','Free','Real',ipRV,3*nsAtom)
           iDo_dDipM=0
           Call GF_on_the_fly(iDo_dDipM)
*
       End If
*                                                                      *
************************************************************************
*                                                                      *
      If (Stop.or.do_printcoords)
     &   Call CollapseOutput(0,'Geometry section')
*     Call QExit('DstInf')
      Return
      End
