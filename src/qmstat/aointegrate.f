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
      Subroutine AOIntegrate(iCStart,nBaseQ,nBaseC,Ax,Ay,Az,nCnC_C
     &          ,iQ_Atoms,nAtomsCC,ipAOint,ipAOintpar,iV2,N,lmax
     &          ,Inside)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "integral.fh"
#include "WrkSpc.fh"
**Taked out by Valera #include "MolcasQmStat.fh"
#include "lenin.fh"

      Dimension V2(MxBasC,MxOrb_C)
      Dimension nCnC_C(MxBasC)
      Dimension Sint(MxBas,MxBasC),SintPar(MxBas,MxBasC),Rot(3,3)
      Dimension Inside(MxAt,3)
      Character Snack*30,BsLbl*4000
      Logical PrEne,PrOcc,Inside

*--------------------------------------------------------------------------*
* Call Transrot. There we compute the rotation matrix for the classical    *
* water under consideration. Used later.                                   *
*--------------------------------------------------------------------------*
      Call TransRot(Cordst,N+1,Rot,Dx,Dy,Dz,Ax,Ay,Az)
      If(iPrint.ge.17) then
        Write(6,*)
        Write(6,*)'ROTATION MATRIX, Molecule ',N/nCent
        Write(6,*)(Rot(1,k),k=1,3)
        Write(6,*)(Rot(2,k),k=1,3)
        Write(6,*)(Rot(3,k),k=1,3)
      Endif
*--------------------------------------------------------------------------*
* Call OrbRot2. Given the rotation matrix (Rot) and the original MO-       *
* coefficients, we transform them to new MO-coefficients. V2 is on input   *
* the original MO-coefficients (stored in V3), and on output the rotated.  *
*--------------------------------------------------------------------------*
      Do 5201, iOrS=1,iOrb(2) !Collect original MO-coeff.
        Do 5202, iBaS=1,nBaseC
          V2(iBaS,iOrS)=V3(iBaS,iOrS)
5202    Continue
5201  Continue
      Call OrbRot2(Rot,V2,iQn,iOrb(2),lMax,nCnC_C)
      kaunt=0
      Do 5211, iMO=1,iOrb(2) !Store the rotated in vector for
        Do 5212, iBa=1,nBaseC !later convinience.
          Work(iV2+kaunt)=V2(iBa,iMO)
          kaunt=kaunt+1
5212    Continue
5211  Continue
      If(iPrint.ge.25) then  !Optional print-out.
        PrOcc=.false.
        PrEne=.false.
        Write(snack,'(A,I3)')'Rotated orbitals for water ',N/ncent
        Call GetMem('PrCMO','Allo','Real',ipPPP,nBaseC*iOrb(2))
        kauntadetta=0
        Do 525, i=1,iOrb(2)
          Do 526, j=1,nBaseC
            Work(ipPPP+kauntadetta)=V2(j,i)
            kauntadetta=kauntadetta+1
526       Continue
525     Continue
        Call NameRun('WRUNFIL')
        Call Get_cArray('Unique Basis Names',BsLbl,LENIN4*nBaseC)
        Call Primo(Snack,PrOcc,PrEne,Dummy,Dummy,1,nBaseC
     &         ,iOrb(2),BsLbl,Dummy,Dummy,Work(ipPPP),3)
        Call GetMem('PrCMO','Free','Real',ipPPP,nBaseC*iOrb(2))
      Endif
      Do 531, m=1,lMax  !New basis function origo definied.
        x=0
        y=0
        z=0
        Do 541, j=1,3
          x=x+Rot(1,j)*SavOri(j,m)
          y=y+Rot(2,j)*SavOri(j,m)
          z=z+Rot(3,j)*SavOri(j,m)
541     Continue
        CasOri(1,m)=x+Dx
        CasOri(2,m)=y+Dy
        CasOri(3,m)=z+Dz
531   Continue
*----------------------------------------------------------------------*
* Compute overlap between the contracted basis functions on the water  *
* molecule presently studied and the QM-molecule.                      *
*----------------------------------------------------------------------*
      Do i=1,nBaseQ
        Do j=1,nBaseC
          Sint(i,j)=0
          SintPar(i,j)=0
        Enddo
      Enddo
      Call ContractOvl(Sint,SintPar,nBaseQ,nBaseC
     &    ,N,nCent,iEl,iQ_Atoms,nAtomsCC,iPrint,Inside)
      !To be able to use the fast matrix multiplication routine DGEMM_,
      !we have to put the Sint (and Sintpar) matrices in vector form.
      !In the future we might 'cut out the middle-man' and already
      !above put the overlap matrix in vector shape.
      kaunt=0
      Do 547, iC=1,nBaseC
        Do 548, iQ=1,nBaseQ
          Work(ipAOint+kaunt)=Sint(iQ,iC)
          kaunt=kaunt+1
548     Continue
547   Continue

      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(iCStart)
        Call Unused_integer(ipAOintpar)
      End If
      End
