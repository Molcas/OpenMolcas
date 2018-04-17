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
* Copyright (C) Anders Ohrn                                            *
************************************************************************
*  OrbRot2
*
*> @brief
*>   Rotate orbitals of the solvent
*> @author A. Ohrn
*>
*> @details
*> Given a rotation matrix \f$ R \f$ it is easy to rotate the p-orbitals since
*> they behave just like the three axes. The d-orbitals are more
*> intricate.
*>
*> To obtain a representation
*> of the d-orbital \f$ \mathrm{d}_{xy} \f$ we perform the outer product \f$ \mathrm{d}_x \mathrm{d}_y + \mathrm{d}_y \mathrm{d}_x \f$,
*> which generates a two-dimensional matrix. Now rotate each px and py
*> and perform this multiplication. We obtain a symmetric matrix from
*> whose elements we can compute how the \f$ \mathrm{d}_{xy} \f$ transforms when
*> rotated. The same is done for the other d-orbitals and we can
*> construct a transformation matrix, which we apply to the MO-coeff.
*> Other ways to look at it are available, but the present one can be
*> seen as a generation of table 3 in \cite Iva1996-JPC-100-6342.
*> The present method is efficient
*> with no trigonometric functions and ample use of BLAS_UTIL.
*>
*> @note
*> ::Qfread as well as ::transrot must precede.
*>
*> @param[in] Rot  The rotation matrix
*> @param[in] Cmo  The MO-coefficients
*> @param[in] iQ   The angular type of the \f$ i \f$ -th basis (observe, not the \f$ i \f$ -th basis function, see givemeinfo.f)
*> @param[in] iOrb Number of orbitals
*> @param[in] lMax Number of bases (not basis functions), see qfread.f
*> @param[in] nCnC Number of contracted basis functions of same type as the \f$ i \f$ -th basis. For example 7s4p will have vector ``7,7,7,7,7,7,7,4,4,4,4``.
************************************************************************
      Subroutine OrbRot2(Rot,Cmo,iQ,iOrb,lMax,nCnC)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "warnings.fh"

      Dimension Rot(3,3),Cmo(MxBasC,MxOrb_C),iQ(MxBasC),nCnC(MxBasC)

      Dimension Px(3),Py(3),Pz(3),Resx(3),Resy(3),Resz(3)
      Dimension Dzz(6),Dxy(6),Dxz(6),Dyz(6),Dxmin(6),Dzer(6)
     &,Dkopia(6),Unit(6)

      Dimension PBlock(3,3),DBlock(5,5)

      Logical NewIq

      Data Px/1.0d0,0.0d0,0.0d0/
      Data Py/0.0d0,1.0d0,0.0d0/
      Data Pz/0.0d0,0.0d0,1.0d0/
      Data Dzer/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
      Data Unit/-1.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,-1.0d0/

*
*-- Dgemv multiply the matrix Rot on the vector Px. For details about
*   various parameter see the source code for the routine; it is very
*   detailed. Here we get the p-orbitals which also are used to obtain
*   the transformation of the other orbitals.
*
      Call dGeMV_('N',ithree,ithree,ONE,Rot,ithree
     &,Px,ione,ZERO,Resx,ione)
      Call dGeMV_('N',ithree,ithree,ONE,Rot,ithree
     &,Py,ione,ZERO,Resy,ione)
      Call dGeMV_('N',ithree,ithree,ONE,Rot,ithree
     &,Pz,ione,ZERO,Resz,ione)

*
*-- Construct the p-block. Easy since they are linear polynomials in
*   the cartesian R-matrix.
*
      PBlock(1,1)=Resx(1)
      PBlock(2,1)=Resx(2)
      PBlock(3,1)=Resx(3)
      PBlock(1,2)=Resy(1)
      PBlock(2,2)=Resy(2)
      PBlock(3,2)=Resy(3)
      PBlock(1,3)=Resz(1)
      PBlock(2,3)=Resz(2)
      PBlock(3,3)=Resz(3)

*
*-- Generation of quadratic polynomial of R-matrix elements for d_xy.
*
      call dcopy_(isix,Dzer,ione,Dkopia,ione)
      Call Dspr2('L',ithree,ONE,Resx,ione,Resy,ione,Dzer)
      call dcopy_(isix,Dzer,ione,Dxy,ione)
      call dcopy_(isix,Dkopia,ione,Dzer,ione)

*
*-- Generation of quadratic polynomial of R-matrix elements for d_xz.
*
      call dcopy_(isix,Dzer,ione,Dkopia,ione)
      Call Dspr2('L',ithree,ONE,Resx,ione,Resz,ione,Dzer)
      call dcopy_(isix,Dzer,ione,Dxz,ione)
      call dcopy_(isix,Dkopia,ione,Dzer,ione)

*
*-- Generation of quadratic polynomial of R-matrix elements for d_yz.
*
      call dcopy_(isix,Dzer,ione,Dkopia,ione)
      Call Dspr2('L',ithree,ONE,Resy,ione,Resz,ione,Dzer)
      call dcopy_(isix,Dzer,ione,Dyz,ione)
      call dcopy_(isix,Dkopia,ione,Dzer,ione)

*
*-- Generation of quadratic polynomial of R-matrix elements for d_xx-yy.
*
      call dcopy_(isix,Dzer,ione,Dkopia,ione)
      Call Dspr('L',ithree,-1.0d0*ONE,Resy,ione,Dzer)
      call dcopy_(isix,Dzer,ione,Dxmin,ione)
      call dcopy_(isix,Dkopia,ione,Dzer,ione)
      Call Dspr('L',ithree,ONE,Resx,ione,Dxmin)
      call dcopy_(isix,Dkopia,ione,Dzer,ione)

*
*-- Generation of quadratic polynomial of R-matrix elements for d_3zz-1.
*
      call dcopy_(isix,Unit,ione,Dkopia,ione)
      Call Dspr('L',ithree,THREE,Resz,ione,Unit)
      call dcopy_(isix,Unit,ione,Dzz,ione)
      call dcopy_(isix,Dkopia,ione,Unit,ione)

*
*-- Construct the d-block. This requires knowledge of how Molcas
*   handles it d-orbitals, especially its normalization, which is the
*   reason for the constants r1,r2,r3,r4.
*
      r1=1.0d0
      r2=sqrt(dble(3))
      r3=1.0d0
      r4=1.0d0/sqrt(dble(3))
      DBlock(1,1)=r1*Dxy(2) !dxy-->dxy
      DBlock(2,1)=r1*Dxy(5) !dxy-->dyz
      DBlock(3,1)=(r1/r4)*0.5*Dxy(6) !dxy-->d3zz-1
      DBlock(4,1)=r1*Dxy(3) !dxy-->dxz
      DBlock(5,1)=(r1/r3)*(Dxy(1)+0.5*Dxy(6)) !dxy-->dxx-yy
      DBlock(1,2)=r1*Dyz(2) !dyz-->dxy
      DBlock(2,2)=r1*Dyz(5) !dyz-->dyz
      DBlock(3,2)=(r1/r4)*0.5*Dyz(6) !dyz-->d3zz-1
      DBlock(4,2)=r1*Dyz(3) !dyz-->dxz
      DBlock(5,2)=(r1/r3)*(Dyz(1)+0.5*Dyz(6)) !dyz-->dxx-yy
      DBlock(1,3)=r4*Dzz(2) !d3zz-->dxy
      DBlock(2,3)=r4*Dzz(5) !d3zz-->dyz
      DBlock(3,3)=r1*0.5*Dzz(6) !d3zz-->d3zz
      DBlock(4,3)=r4*Dzz(3) !d3zz-->dxz
      DBlock(5,3)=(r1/r2)*(Dzz(1)+0.5*Dzz(6)) !d3zz-->dxx-yy
      DBlock(1,4)=r1*Dxz(2) !dxz-->dxy
      DBlock(2,4)=r1*Dxz(5) !dxz-->dyz
      DBlock(3,4)=(r1/r4)*0.5*Dxz(6) !dxz-->d3zz-1
      DBlock(4,4)=r1*Dxz(3) !dxz-->dxz
      DBlock(5,4)=(r1/r3)*(Dxz(1)+0.5*Dxz(6)) !dxz-->dxx-yy
      DBlock(1,5)=r3*Dxmin(2) !dxx-yy-->dxy
      DBlock(2,5)=r3*Dxmin(5) !dxx-yy-->dyz
      DBlock(3,5)=r2*0.5*Dxmin(6) !dxx-yy-->d3zz
      DBlock(4,5)=r3*Dxmin(3) !dxx-yy-->dxz
      DBlock(5,5)=r1*(Dxmin(1)+0.5*Dxmin(6)) !dxx-yy-->dxx-yy

*
*-- With the proper number of blocks at hand, we make transformation.
*   The brute force way is to construct the entire transformation matrix
*   and make a matrix multiplication. Since so many elements in that
*   matrix are zero, more effecient ways to transform are available.
*   That is what is used below. The formulas follow from considering the
*   transformation matrix multiplied with the CMO-matrix. Multio
*   importante is it to observe how the basis functions in Molcas are
*   ordered!
*
      Do 10, i=1,iOrb
        In=1   !OBSERVE!!! WE ARE ASSUMING THAT THE FIRST BASIS IS OF
               !S-TYPE!!! IF YOU ARE IMPLEMENTING SYMMETRY, THIS
        Do 20, j=2,lMax  !MIGHT NOT BE A VALID ASSUMPTION SO THEN THE
          In=In+1       !NEWIQ-CONSTRUCT BELOW MUST BE ALTERED!!!
          IqSuckOut=iQ(j)
          NewIq=iQ(j).ne.iQ(j-1)
          If(Newiq) then !This if-clause controls the jumping when new
            In=In+(2*iQ(j-1)-2)*nCnC(j-1) !angular basis function
          Endif                           !appears.
          If(IqSuckOut.eq.1) then   !This is s-function
            Go To 20
          ElseIf(IqSuckOut.eq.2) then  !This is p-function
            iSkutt=nCnC(j)
            Ctemp1=PBlock(1,1)*Cmo(In,i)+PBlock(1,2)*Cmo(In+iSkutt,i)
     &           +PBlock(1,3)*Cmo(In+2*iSkutt,i)
            Ctemp2=PBlock(2,1)*Cmo(In,i)+PBlock(2,2)*Cmo(In+iSkutt,i)
     &           +PBlock(2,3)*Cmo(In+2*iSkutt,i)
            Ctemp3=PBlock(3,1)*Cmo(In,i)+PBlock(3,2)*Cmo(In+iSkutt,i)
     &           +PBlock(3,3)*Cmo(In+2*iSkutt,i)
            Cmo(In,i)=Ctemp1
            Cmo(In+iSkutt,i)=Ctemp2
            Cmo(In+2*iSkutt,i)=Ctemp3
          ElseIf(IqSuckOut.eq.3) then  !This is d-function
            iSkutt=nCnC(j)
            Ctemp1=DBlock(1,1)*Cmo(In,i)+DBlock(1,2)*Cmo(In+iSkutt,i)
     &            +DBlock(1,3)*Cmo(In+2*iSkutt,i)
     &            +DBlock(1,4)*Cmo(In+3*iSkutt,i)
     &            +DBlock(1,5)*Cmo(In+4*iSkutt,i)
            Ctemp2=DBlock(2,1)*Cmo(In,i)+DBlock(2,2)*Cmo(In+iSkutt,i)
     &            +DBlock(2,3)*Cmo(In+2*iSkutt,i)
     &            +DBlock(2,4)*Cmo(In+3*iSkutt,i)
     &            +DBlock(2,5)*Cmo(In+4*iSkutt,i)
            Ctemp3=DBlock(3,1)*Cmo(In,i)+DBlock(3,2)*Cmo(In+iSkutt,i)
     &            +DBlock(3,3)*Cmo(In+2*iSkutt,i)
     &            +DBlock(3,4)*Cmo(In+3*iSkutt,i)
     &            +DBlock(3,5)*Cmo(In+4*iSkutt,i)
            Ctemp4=DBlock(4,1)*Cmo(In,i)+DBlock(4,2)*Cmo(In+iSkutt,i)
     &            +DBlock(4,3)*Cmo(In+2*iSkutt,i)
     &            +DBlock(4,4)*Cmo(In+3*iSkutt,i)
     &            +DBlock(4,5)*Cmo(In+4*iSkutt,i)
            Ctemp5=DBlock(5,1)*Cmo(In,i)+DBlock(5,2)*Cmo(In+iSkutt,i)
     &            +DBlock(5,3)*Cmo(In+2*iSkutt,i)
     &            +DBlock(5,4)*Cmo(In+3*iSkutt,i)
     &            +DBlock(5,5)*Cmo(In+4*iSkutt,i)
            Cmo(In,i)=Ctemp1
            Cmo(In+iSkutt,i)=Ctemp2
            Cmo(In+2*iSkutt,i)=Ctemp3
            Cmo(In+3*iSkutt,i)=Ctemp4
            Cmo(In+4*iSkutt,i)=Ctemp5
          Else  !Here we go if non-implemented angular quantum number
                !appears.
            Write(6,*)
            Write(6,*)' ERROR in OrbRot2. I''m not ready for f-orbitals'
            Call Quit(_RC_GENERAL_ERROR_)
          Endif
20      Continue
10    Continue

      Return
      End
