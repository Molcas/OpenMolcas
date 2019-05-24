************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1999, Coen de Graaf                                    *
*               1999, Anders Bernhardsson                              *
*               1999, Roland Lindh                                     *
*               2018, Jesper Norell                                    *
************************************************************************
      Subroutine Molden_DysOrb(iUHF,filename,ENE,OCC,CMO,NDO,NZ)
************************************************************************
*                                                                      *
*     Object: to generate MOLDEN files for Dyson orbitals              *
*                                                                      *
*     Modified from molden_interface by Jesper Norell 2018             *
*     (Since requirements for e.g. number of orbitals and              *
*     treament of symmetry differs from "normal" MOs)                  *
*                                                                      *
*     Authors: Coen de Graaf, Anders Bernardsson and R. Lindh, 1999    *
*                                                                      *
************************************************************************
      use Real_Spherical
      implicit real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
*
*
c      Parameter (MaxOrb_Molden=400, MaxOrb_Do=100)
      Parameter (EorbThr = 50.D0 )
      Real*8 Coor(3,mxdc),Znuc(mxdc)
      Character shelllabel(7)
      Character*(LENIN) AtomLabel(mxdc)
      Character*(LENIN8) label(MaxBfn+MaxBfn_Aux)
      Character*8 MO_Label(maxbfn)
      Parameter (nNumber=61)
      Character Number(nNumber)
      Integer ibas_lab(mxdc), nOrb(8)
      Character*(LENIN8+1) gtolabel(maxbfn)
      Real*8 r_Norm(maxbfn)
      Character*(*) Filename
      Character Env*8
      Logical Exist,y_cart,y_sphere, Found, Reduce_Prt
      External Reduce_Prt
      Character*100 Supername,Get_SuperName
      External Get_SuperName
      data shelllabel /'s','p','d','f','g','h','i'/
      data number /'1','2','3','4','5','6','7','8','9','0',
     &             'a','b','c','d','e','f','g','h','i','j',
     &             'k','l','m','n','o','p','q','r','s','t',
     &             'u','v','w','x','y','z','A','B','C','D',
     &             'E','F','G','H','I','J','K','L','M','N',
     &             'O','P','Q','R','S','T','V','W','X','Y',
     &             'Z'/
      data iRc/0/
      save iRc

      ! Jesper
      Integer NDO,NZ
      Dimension ENE(NDO)
      Dimension OCC(NDO)
      Dimension CMO(NZ*NDO)
      Dimension DESYM(NZ,NZ)
*
*     Statement function
*
      CC(ix,iy,iz)=SQRT(DblFac(2*ix-1)*DblFac(2*iy-1)*DblFac(2*iz-1))
*
      if(iRc.eq.1) Return
*
* Do nothing within numerical_gradient
      SuperName=Get_Supername()
      If (SuperName == 'numerical_gradient') Then
        Return
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Get the print level
*
      iPL=iPrintLevel(-1)
      jPL=iPL
      If (Reduce_Prt().and.iPL.lt.3) jPL=0
*                                                                      *
************************************************************************
*                                                                      *
      If (MolWgh.eq.1) Then
         If (jPL.ge.2) Then
            Write (6,*) 'Molden_Interface: '
     &                //'Unsupported normalization,Molwgh=1!'
         End If
         iRc=1
         Return
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call GetEnvf('MOLCAS_MOLDEN',Env)
c      If (Env.eq.' '.or.Env.eq.'OFF') Then
      If (Env.eq.'OFF') Then
         If (jPL.ge.2) Then
            Write(6,*)
            Write(6,*) ' Molden files will not be produced'
            Write(6,*)
         End If
         iRC=1
         Return
      End If
cVV: current version of Molden has no clear limit for MaxOrb
c      If (MaxOrb.gt.MaxOrb_Molden) Then
c         If (jPL.ge.2) Then
c            Write(6,*)
c            Write(6,*) ' Molden_Interface: W A R N I N G !!!!'
c            Write(6,*)
c            Write(6,*) ' No Molden input file will be generated!'
c            Write(6,*)
c            Write(6,*) ' Calculation exceeds the max number of orbitals'
c     &               //' allowed for MOLDEN. To change this modify the '
c            Write(6,*) ' parameter MaxOrb_Molden in  src/util/molden_in'
c     &               //'terface.f and follow the instructions in Molden'
c            Write(6,*) ' on how to modify the parameter MaxOrb.'
c            Write(6,*)
c         End If
c         iRC=1
c         Return
c      End If
*                                                                      *
************************************************************************
*                                                                      *
      Check_CMO=Zero
      Check_Energy=Zero
      Check_Occupation=Zero
*                                                                      *
************************************************************************
*                                                                      *
      Call f_Inquire('RUNFILE',Exist)
      If (.Not.Exist) then
      iRC=1
      Return
      Endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Read the characteristics of all different basis sets,
*     provide each atom with a nuclear charge and establish
*     a link between an atom and its basis set ---
*
*     NOTICE!!!
*     This call will also fill info.fh and the dynamic storage in
*     Work(ipInf)
*
      Call Inter1       (AtomLabel,iBas_Lab,Coor,Znuc,nAtom,ipInf)
      Call Qpg_iArray('nOrb',Found,nData)
      If (Found) Then
         Call Get_iArray('nOrb',nOrb,nData)
      Else
         Call iCopy(nIrrep,nBas,1,nOrb,1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      iAngMx_Valence=0
      Do iCnttp = 1, nCnttp
         If (.Not.AuxCnttp(iCnttp) .and.
     &       .Not.FragCnttp(iCnttp) ) Then
            nTest=nVal_Shells(iCnttp)-1
            iAngMx_Valence=Max(iAngMx_Valence,nTest)
         End If
      End Do
      If (iAngMx_Valence.gt.4) Then
        If (jPL.ge.2) Then
           Write(6,*) 'Sorry, Molden does not know how to handle'
           Write(6,*) 'functions with angular momentum larger than g'
        End If
        Go To 999
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Unnormalize contraction coefficients for the valence shells
*
      Do iCnttp=1,nCnttp
        If (.Not.(AuxCnttp(iCnttp).or.FragCnttp(iCnttp))) Then
         Do l=0,nVal_Shells(iCnttp)-1
          ishell=ipVal(iCnttp)+l
          Call Unnrmlz(Work(ipExp(ishell)),nexp(ishell),
     &                 Work(ipCff(ishell)),nbasis(ishell),l)
         End Do
        End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Compute memory requirements and allocate memory
*
      nB=0
      Do iS=0,nirrep-1
       nB=nB+nBas(is)
      End Do
      Call GetMem('ICENT','ALLO','INTE',ipCent,8*nB)
      Call GetMem('IPHASE','ALLO','INTE',ipPhase,8*nB)
      Call GetMem('nCENT','ALLO','INTE',ipCent2,nB)
      Call GetMem('ICENTER','ALLO','INTE',ipCent3,nB)
      Call GetMem('CMO2','ALLO','REAL',ipC2,nB**2)
      Call GetMem('VECTOR','ALLO','REAL',ipV,nB**2)
      call dcopy_(nB**2,[Zero],0,Work(ipV),1)
      If (iUHF.eq.1) Then
         Call GetMem('CMO2','ALLO','REAL',ipC2_ab,nB**2)
         Call GetMem('VECTOR','ALLO','REAL',ipV_ab,nB**2)
         call dcopy_(nB**2,[Zero],0,Work(ipV_ab),1)
      Else
         ipC2_ab=ip_Dummy
         ipV_ab =ip_Dummy
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Open input file for MOLDEN
*
      MF=9
      call molcas_open(MF,filename)
*                                                                      *
************************************************************************
*                                                                      *
*     Write the atom information in MOLDEN format to unit MF
*
*
      y_cart=.false.
      y_sphere=.false.
      Do iCnttp=1,nCnttp
        If (AuxCnttp(iCnttp).or.FragCnttp(iCnttp)) Go To 995
        Do iCntr=1,nCntr(iCnttp)
          Do l=0,nVal_Shells(iCnttp)-1
*           Test for the appearance of cartesian functions with l=2,3,4
            ishell=ipVal(iCnttp)+l
            if ((l.ge.2).and.(.not.y_cart)) Then
              if (.not.transf(ishell)) y_cart=.true.
            End If
            if ((l.ge.2).and.(.not.y_sphere)) Then
              if (transf(ishell)) y_sphere=.true.
            end if
            if (y_sphere.and.y_cart) Then
              If (jPL.ge.2) Then
                 Write (6,*)
                 Write (6,*) 'Failed to generate input file to MOLDEN'
                 Write (6,*) 'No mixing allowed of spherical and',
     &                       ' cartesian d, f, g-functions'
              End If
              Go to 991
            End If
          End Do
        End Do
 995    Continue
      End Do
      Write (MF,'(A)') '[MOLDEN FORMAT]'
*                                                                      *
************************************************************************
*                                                                      *
*     Write atomic information
*
      Write (MF,'(A)') '[N_ATOMS]'
      Write (MF,*) natom
      Write (MF,'(A)') '[ATOMS] (AU)'
      Do iatom=1,natom
        Write (MF,99) AtomLabel(iatom),iatom,int(Znuc(iatom)),
     &                (coor(i,iatom),i=1,3)
 99     Format(A,2(3x,I4),5x,3F16.8,3I4)
      End Do
      If (.not.y_cart) Then
        If (iAngMx.gt.1) Write (MF,'(A)') '[5D]'
        If (iAngMx.gt.2) Write (MF,'(A)') '[7F]'
        If (iAngMx.gt.3) Write (MF,'(A)') '[9G]'
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Write Gaussian basis set information to MOLDEN input file.
*
      Write (MF,'(A)') '[GTO] (AU)'
*
*     Read exponents and contraction coefficients of each unique basis.
*     Write the present basis set (iCnttp) to the molden.input file for
*     the appropriate atoms.
*     Moreover, a list is constructed which contains a label for each
*     GTO (gtolabel). This list follows the MOLDEN order of GTOs.
*     Later this list will be used in the transformation of sabf (the
*     symmetry adapted basis functions).
*
      iatom=0
      mdc=0
      kk=0
*
      Do iCnttp=1,nCnttp             ! loop over unique basis sets
        If (AuxCnttp(iCnttp).or.FragCnttp(iCnttp)) Go To 996
*
        Do iCntr=1,nCntr(iCnttp)     ! loop over sym. unique centers
          mdc=mdc+1
          nDeg=nIrrep/nStab(mdc)
          Do iDeg=1,nDeg             ! loop over centers
            iAtom=iAtom+1
            Write (MF,'(I4)') iAtom
*
            Do l=0,nVal_Shells(iCnttp)-1
              ishell=ipVal(iCnttp)+l
              If (nBasis(iShell).gt.nNumber) Then
                 Write (6,*) 'Interf: too many contracted functions!'
                 Write (6,*) 'nBasis(iShell)=',nBasis(iShell)
                 Call Abend()
              End If
*
*             Iterate over each contracted GTO
*
              Do icontr=1,nBasis(ishell)
*
*               Find the number of exponents with non-zero exponents
*
                isegm=0
                Do iprim=1,nExp(ishell)
                  coeff=
     &             Work(ipCff(ishell)+(icontr-1)*nExp(ishell)+iprim-1)
                  If (coeff.ne.Zero) Then
                    isegm=isegm+1
                  End If
                End Do
*
                Write (MF,'(3x,A1,I4)') shelllabel(l+1),isegm
*
*               Write exponents and contraction coefficients.
*
                Do iprim=1,nExp(ishell)
                  coeff=
     &             Work(ipCff(ishell)+(icontr-1)*nExp(ishell)+iprim-1)
                  prim=work(ipExp(ishell)+iprim-1)
                  If (coeff.ne.Zero) Then
                    Write (MF,'(E17.9,E17.9)') prim,coeff
                  End If
                End Do
*
*               Construction of gtolabel
*               Molden order: for p-functions: px(1), py(1), pz(1),
*                                              px(2), py(2), pz(2), etc.
*               for d-, and f-functions: similar
*
                If (l.eq.0) Then
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'01s     '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                End If
                If (l.eq.1) Then
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'02px    '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'02py    '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'02pz    '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                End If
                If ((l.eq.2).and.(.not.y_cart)) Then
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'03d00   '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'03d01+  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'03d01-  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'03d02+  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'03d02-  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                End If
                If ((l.eq.2).and.(y_cart)) Then
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'d020000 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(2,0,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'d000200 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,2,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'d000002 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,0,2)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'d010100 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(1,1,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'d010001 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(1,0,1)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'d000101 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,1,1)
                  iWork(ipCent3+kk-1)=iAtom
                End If
                If ((l.eq.3).and.(.not.y_cart)) Then
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'04f00   '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'04f01+  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'04f01-  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'04f02+  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'04f02-  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'04f03+  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'04f03-  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                End If
                If ((l.eq.3).and.(y_cart)) Then
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f030000 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(3,0,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f000300 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,3,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f000003 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,0,3)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f010200 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(1,2,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f020100 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(2,1,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f020001 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(2,0,1)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f010002 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(1,0,2)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f000102 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,1,2)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f000201 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,2,1)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'f010101 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(1,1,1)
                  iWork(ipCent3+kk-1)=iAtom
                End If
                If ((l.eq.4).and.(.not.y_cart)) Then
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'05g00   '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'05g01+  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'05g01-  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'05g02+  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                 iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'05g02-  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'05g03+  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'05g03-  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'05g04+  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'05g04-  '//
     &                         number(icontr)
                  r_Norm(kk)=1.0D0
                  iWork(ipCent3+kk-1)=iAtom
                End If
                If ((l.eq.4).and.(y_cart)) Then
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g040000 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(4,0,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g000400 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,4,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g000004 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,0,4)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g030100 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(3,1,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g030001 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(3,0,1)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g010300 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(1,3,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g000301 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,3,1)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g010003 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(1,0,3)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g000103 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,1,3)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g020200 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(2,2,0)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g020002 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(2,0,2)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g000202 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(0,2,2)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g020101 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(2,1,1)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g010201 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(1,2,1)
                  iWork(ipCent3+kk-1)=iAtom
                  kk=kk+1
                  gtolabel(kk)=AtomLabel(iAtom)//'g010102 '//
     &                         number(icontr)
                  r_Norm(kk)=CC(1,1,2)
                  iWork(ipCent3+kk-1)=iAtom
                End If
              End Do
            End Do
            Write (MF,'(A)') ' '
          End Do
        End Do
 996    Continue
      End Do
      kk_Max=kk
      If (nB.gt.kk_max) Then
         If (jPL.ge.2) Then
            Write(6,*) 'Molden_Interface: nB.gt.kk_max'
            Write(6,*) 'nB,kk_Max=',nB,kk_Max
         End If
         Go To 998
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
      nTot=0
      nTot2=0
      Do iS=0,nIrrep-1
         nTot=nTot+nBas(iS)
         nTot2=nTot2+nBas(iS)**2
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*      Back 'transformation' of the symmetry adapted basis functions.
*      Probably somewhat clumsy, but it seems to work.If someone
*      knows a more elegant way to do it, please improve this part!
*
*      PART 1: Obtain symmetry information (soout), construct a label
*              for each sabf, which will be used in part 2 to find the
*              corresponding GTO in the MOLDEN list by comparing with
*              gtolabel
*
*      nB       --- Total number of contracted basis functions
*      ipcent2  --- degeneracy of a basis function
*      ipCent   --- centres over which the basis function is
*                   delocalized
*      ipPhase  --- phase of the AO in the linear combination
*
      Call icopy(8*nB,[0],0,iWork(ipPhase),1)
      Call icopy(8*nB,[0],0,iWork(ipCent),1)
      Call SOout(label,iWork(ipCent),iWork(ipPhase))
      ipc=0
      Do iContr=1,nB
        iWork(ipCent2+iContr-1)=0
        Do k=1,8
          If (iWork(ipCent+ipc).ne.0)
     &       iWork(ipcent2+iContr-1)=iWork(ipCent2+iContr-1)+1
          ipc=ipc+1
        End Do
      End Do
*

***********************************************************************
*
*      Construct a matrix which inverses the symmetry transformation
*      of the basis fuctions (DESYM)
*
      IF (nIrrep.LT.2) THEN
          GOTO 500
      END IF

      DESYM=0.0D0
      iBtot=0
      iCount=0
      Do iIrrep=0,nIrrep-1 ! For all the irreps of symmetrized functions

        Do iB=1,nBas(iIrrep) ! For each symmetrized function in the
          iBtot=iBtot+1      ! irrep
          iSymcent=0

          If (iB.eq.1) Then ! Counter within functions with the
            iCount=1        ! same label (e.g. H1 1s 1, H1 1s 2, ...)
          else
            If (label(iBtot-1).eq.label(iBtot)) Then
              iCount=iCount+1
            else
              iCount=1
            End If
          End If

          Write (MO_Label(iBtot),'(A3)') lirrep(iIrrep) ! Save label
          ! for later printing

          Do iGTO=1,nB ! Find matching Gaussian (GTO)

            If (gtolabel(iGTO).eq.label(iBtot)//number(iCount)) Then
                ! GTO Found
                idx1=8*(iBtot-1)+iSymcent
                idx2=(iBtot-1)
                coeff=DBLE(iWork(ipPhase+idx1))*r_Norm(iGTO)
     &                / SQRT(DBLE(iWork(ipCent2+idx2)))
                DESYM(iBtot,iGTO)=DESYM(iBtot,iGTO)+coeff
                iSymcent=iSymcent+1 ! Keep track of symmetry centre
            End If ! GTO finding

          End Do ! iGTO=1,nB
        End Do ! iB=1,nBas(iIrrep)
      End Do ! iIrrep=0,nIrrep-1

************************************************************************
*                                                                      *
*
*                                                                      *
************************************************************************
*                                                                      *
*      Dump vector in the molden.input file
*

 500  Continue
      Write (MF,'(A)') '[MO]'

      ! For all Dyson orbitals
      DO I=1,NDO

      ! No symmetry = no backtransformation of basis functions needed
       IF (nIrrep.LT.2) THEN
        Write (MF,'(A,I0,A)') 'Sym= ',I,"a"
        Write (MF,103) ENE(I)
        Write (MF,'(A)') 'Spin= Alpha'
        Write (MF,104) OCC(I)
        DO J=1,NB
         WRITE(MF,100) J,CMO((I-1)*NB+J)
        END DO

      ! Symmetry = perform a backtransformation from symmetrized to
      ! original basis functions
       ELSE
        ! Find the correct symmetry label from the symmetrized functions
        idx=(I-1)*NB+1
        iLabel=MAXLOC(ABS(CMO(idx:idx+nB-1)),1)
        Write (MF,'(A,I0,A)') 'Sym= ',I,MO_Label(iLabel)
        Write (MF,103) ENE(I)
        Write (MF,'(A)') 'Spin= Alpha'
        Write (MF,104) OCC(I)

        ! Now backtransform the orbitals
        DO iGTO=1,NB ! For each molden GTO
         COEF=0
         DO iB=1,NB ! Gather contributions from all symmetrized
                    ! basis functions
          COEF=COEF+DESYM(iB,iGTO)*CMO((I-1)*NB+iB)
         END DO
         WRITE(MF,100) iGTO,COEF
        END DO
        CONTINUE

       END IF ! End of symmetry backtransformation

      END DO ! I=1,NDO (i.e. Dyson orbitals)



*                                                                      *
************************************************************************
*                                                                      *
 991  Continue
      Call GetMem('ICENT','FREE','INTE',ipCent,8*nB)
      Call GetMem('IPHASE','FREE','INTE',ipPhase,8*nB)
      Call GetMem('nCENT','FREE','INTE',ipCent2,nB)
      Call GetMem('ICENTER','FREE','INTE',ipCent3,nB)
      If (iUHF.eq.1) Then
         Call GetMem('CMO2','FREE','REAL',ipC2_ab,nB**2)
         Call GetMem('VECTOR','FREE','REAL',ipV_ab,nB**2)
      End If
      Call GetMem('CMO2','FREE','REAL',ipC2,nB**2)
      Call GetMem('VECTOR','FREE','REAL',ipV,nB**2)
 998  Close (MF)
 999  Call ClsSew
*
C -------------FORMATS-------------
 100  format(I4,3x,F16.8)
 103  format('Ene= ',F10.4)
 104  format('Occup= ',F10.5)

      Return
      End
