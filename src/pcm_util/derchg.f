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
      Subroutine DerChg(nAt,nAt3,nTs,nS,Eps,IAtm,AtmC,AtmChg,
     &  DM,DerMat,Tessera,Q,Qtot,QDer,
     $  DerTes,DerPunt,DerCentr,DerRad,Der1,Der2,VDer,Sphere,ISphe)
      Implicit Real*8(A-H,O-Z)
      Integer IAtm(nAt)
      Dimension AtmC(3,nAt),AtmChg(nAt)
      Dimension Q(2,*),ISphe(*),QTot(*),QDer(3,nAt,*),Der1(*),Der2(*)
      Dimension Tessera(4,*),Sphere(4,*),VDer(nTs,*)
      Dimension DM(nTs,*),DerMat(nTs,*)
      Dimension DerTes(nTs,NAt,3),DerPunt(nTs,NAt,3,3)
      Dimension DerRad(nS,NAt,3),DerCentr(nS,NAt,3,3)
      Data Zero,One,Two,Four/0.0d0,1.0d0,2.0d0,4.0d0/,XX/0.0d0/,JJ/0/
*
      PI  = Four*ATan(One)
      FPI = Four*PI
      Diag = - 1.0694d0 * Sqrt(FPI) / Two
      Sc_Cond = (Eps - One) / Eps
*
*---- Total charges
*
      Do iTs = 1, nTs
        Qtot(iTs) = Q(1,iTs) + Q(2,iTs)
      endDo
*
*---- Compute the derivative of PCM matrix
*
c     Loop on atoms and coordinates
      Do 100 iAt = 1, nAt
        Do 101 iCoord = 1, 3
          Index = 3 * (iAt-1) + iCoord
c         Conductor-like case
          Call DMat_CPCM(iAt,iCoord,Eps,nTs,nS,nAt,Diag,Tessera,
     &                   DerMat,DerTes,DerPunt,DerCentr,iSphe)

*
*---- First product: DerDM * Qtot
*
          Call PrMatVec(.False.,.False.,DerMat,-One,nTs,nTs,Qtot,Der1)
cpcm_solvent
c      if(iat.eq.5.and.icoord.eq.1) then
c       write(6,'(''In DerChg first contribution for 5, 1'')')
c       do its = 1, nts
c        write(6,'(i4,f20.12)') its,der1(its)
c       enddo
c      endif
cpcm_solvent end
*
*---- Total deriv. of the potential summed up the quantity alread
*---- computed (DerDM*Qtot)
*
          Do 200 iTs = 1, nTs
            Der1(iTs) = Der1(iTs) + Sc_Cond * VDer(iTs,Index)
  200     Continue
*
*---- Last product: - DM^-1 * (V^x + DM^x*q)
*
cpcm_solvent
c      if(iat.eq.5.and.icoord.eq.1) then
c       write(6,'(''In DerChg second contribution for 5, 1'')')
c       do its = 1, nts
c        write(6,'(i4,f20.12)') its,der1(its)
c       enddo
c      endif
cpcm_solvent end
          Call PrMatVec(.False.,.False.,DM,-1.d0,nTs,nTs,Der1,Der2)
          Call FillQDer(nAt,nTs,iAt,iCoord,Der2,QDer)
  101   Continue
  100 Continue
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(nAt3)
        Call Unused_integer_array(IAtm)
        Call Unused_real_array(AtmC)
        Call Unused_real_array(AtmChg)
        Call Unused_real_array(DerRad)
        Call Unused_real_array(Sphere)
      End If
      End
c----------------------------------------------------------------------
      Subroutine DMat_CPCM(iAt,iC,Eps,nTs,nS,nAt,fact,Tessera,
     &           DerMat,DerTes,DerPunt,DerCentr,iSphe)
      Implicit Real*8 (A-H,O-Z)
      Dimension ISphe(*),Tessera(4,*),DerMat(nTs,*)
      Dimension DerTes(nTs,NAt,3),DerPunt(nTs,NAt,3,3)
      Dimension DerCentr(nS,NAt,3,3)
      Data One/1.0d0/,Two/2.0d0/,Four/4.0d0/
C
C     Compute the derivative of the CPCM matrix wrt atom iat, coord. ic
C
C     Loop on tesserae
      Do 10 ITs = 1, NTs
        L = ISPHE(ITs)
        Do 11 JTS = 1, NTs
          LJ = ISPHE(JTS)
C         Diagonal elements
          If(ITs.eq.JTS) then
            DerMat(ITs,ITs) = fact * DERTES(ITs,iAt,IC) /
     &      ( Tessera(4,ITs) * Sqrt(Tessera(4,ITs)) )
          Else
C           Off diagonal elements
            XIJ = Tessera(1,ITs) - Tessera(1,JTS)
            YIJ = Tessera(2,ITs) - Tessera(2,JTS)
            ZIJ = Tessera(3,ITs) - Tessera(3,JTS)
            DIJ = Sqrt( XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ )
            DXIJ=DERPUNT(ITs,iAt,IC,1)+DERCENTR(L,iAt,IC,1)
     &          -DERPUNT(JTS,iAt,IC,1)-DERCENTR(LJ,iAt,IC,1)
            DYIJ=DERPUNT(ITs,iAt,IC,2)+DERCENTR(L,iAt,IC,2)
     &          -DERPUNT(JTS,iAt,IC,2)-DERCENTR(LJ,iAt,IC,2)
            DZIJ=DERPUNT(ITs,iAt,IC,3)+DERCENTR(L,iAt,IC,3)
     &          -DERPUNT(JTS,iAt,IC,3)-DERCENTR(LJ,iAt,IC,3)
            PROD = (XIJ*DXIJ + YIJ*DYIJ + ZIJ*DZIJ) / DIJ**3
            DerMat(ITs,JTS) = - PROD
          EndIf
   11   Continue
   10 Continue
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(Eps)
      End
*------------------------------------------------------------------------
      Subroutine PrMatVec(Dag,DoSym,Mat,f,n,m,Vec,Res)
      Implicit Real*8 (A-H,O-Z)
      Real*8 Mat(n,*),Vec(*),Res(*)
      Logical Dag,DoSym
      Data Zero,Two /0.0d0,2.0d0/
*
*---- Do the matrix vector product: f*Mat(n,m)*Vec(m,1)=Res(n,1)
*---- possibly transposed (if Dag): f*Vec(1,m)*Mat(m,n)=Res(1,n)
*---- If DoSym symmetrize the matrix elements
*
      Do 100 i = 1, n
        Res(i) = Zero
        Do 101 j = 1, m
          If(DoSym) then
            ElM = (Mat(i,j) + Mat(j,i)) / Two
          Else
            If(Dag) ElM = Mat(j,i)
            If(.not.Dag) Elm = Mat(i,j)
          EndIF
          Res(i) = Res(i) + f * ElM * Vec(j)
  101   Continue
  100 Continue
      Return
      End
*------------------------------------------------------------------------
      Subroutine FillQDer(nAt,nTs,iAt,iC,Der,QDer)
      Implicit Real*8 (A-H,O-Z)
      Dimension Der(*),QDer(3,nAt,*)
      Do 100 iTs = 1, nTs
        QDer(iC,iAt,iTs) = Der(iTs)
  100 Continue
      Return
      End
*------------------------------------------------------------------------
      subroutine testq(nAt,nTs,VDer,q,qtot)
      Implicit Real*8 (A-H,O-Z)
      Dimension VDer(nTs,*),Q(2,*),QTot(*)
      Integer Lu
      Lu=1
      Call Molcas_open(Lu,'DerPt.dat')
*     open(1,file='DerPot.dat',status='old',form='formatted')
      do 1132 iAt = 1, nAt
        do 1133 iCoord = 1, 3
          Index = 3 * (iAt-1) + iCoord
          do 1134 iTs = 1, nTs
            read(1,*)VDer(iTs,Index)
 1134     continue
 1133   continue
 1132 continue
      close(1)
      do 1135 iAt = 1, nAt
        do 1136 iCoord = 1, 3
          Index = 3 * (iAt-1) + iCoord
          sum = 0.d0
          do its = 1, nts
            qtot(its) = q(1,its) + q(2,its)
            sum = sum + qtot(its) * VDer(its,index)
          enddo
          write(6,'(''Charges times VDer'',i4,f20.12)') index,sum
 1136   continue
 1135 continue
      return
      end
