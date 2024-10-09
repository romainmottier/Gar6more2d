!! Copyright INRIA. Contributors : Julien DIAZ and Abdelaaziz EZZIANI
!! 
!! Julien.Diaz@inria.fr and Abdelaaziz.Ezziani@univ-pau.fr
!! 
!! This software is a computer program whose purpose is to
!! compute the analytical solution of problems of waves propagation in two 
!! layered media such as
!! - acoustic/acoustic
!! - acoustic/elastodynamic
!! - acoustic/porous
!! - porous/porous,
!! based on the Cagniard-de Hoop method.
!! 
!! This software is governed by the CeCILL license under French law and
!! abiding by the rules of distribution of free software.  You can  use, 
!! modify and/ or redistribute the software under the terms of the CeCILL
!! license as circulated by CEA, CNRS and INRIA at the following URL
!! "http://www.cecill.info". 
!! 
!! As a counterpart to the access to the source code and  rights to copy,
!! modify and redistribute granted by the license, users are provided only
!! with a limited warranty  and the software's author,  the holder of the
!! economic rights,  and the successive licensors  have only  limited
!! liability. 
!! 
!! In this respect, the user's attention is drawn to the risks associated
!! with loading,  using,  modifying and/or developing or reproducing the
!! software by the user in light of its specific status of free software,
!! that may mean  that it is complicated to manipulate,  and  that  also
!! therefore means  that it is reserved for developers  and  experienced
!! professionals having in-depth computer knowledge. Users are therefore
!! encouraged to load and test the software's suitability as regards their
!! requirements in conditions enabling the security of their systems and/or 
!! data to be ensured and,  more generally, to use and operate it in the 
!! same conditions as regards security. 
!! 
!! The fact that you are presently reading this means that you have had
!! knowledge of the CeCILL license and that you accept its terms.
!! ========================================================================

function calccoeff_ap(pp)
  Use m_const
  Use m_phys
implicit none
complex*16,intent(in) :: pp
integer ::Info,IPIV(4)
complex*16 :: k1,kf2,ks2,kpsi2
complex*16 :: Mat(4,4),coeff(4),calccoeff_ap(4)
real*8 :: X(2)
if(real(pp).ne.0.d0) then
   k1=sqrt(pp**2+1/V1**2)
   kf2=sqrt(pp**2+1/Vf2**2)
   ks2=sqrt(pp**2+1/Vs2**2)
   kpsi2=sqrt(pp**2+1/Vpsi2**2)
else
   if (aimag(pp)**2>1/V1**2) then
      k1=ima*sqrt(-pp**2-1/V1**2)
    else
       k1=sqrt(pp**2+1/V1**2)
    end if
    if (aimag(pp)**2>1/Vf2**2) then
       kf2=ima*sqrt(-pp**2-1/Vf2**2)
    else
       kf2=sqrt(pp**2+1/Vf2**2)
    end if
    if (aimag(pp)**2>1/Vs2**2) then
       ks2=ima*sqrt(-pp**2-1/Vs2**2)
    else
       ks2=sqrt(pp**2+1/Vs2**2)
    end if
    if (aimag(pp)**2>1/Vpsi2**2) then
       kpsi2=ima*sqrt(-pp**2-1/Vpsi2**2)
    else
       kpsi2=sqrt(pp**2+1/Vpsi2**2)
    end if
 end if
if (open.eq.1) then
!!!! Open pore transmission conditions
  MAT(1,1)=-k1/rho1
  MAT(2,1)=1
  MAT(3,1)=0
  MAT(4,1)=1

  X=P2(:,1)
  MAT(1,2)=kf2*X(2)+kf2*X(1);
  MAT(2,2)=M2/Vf2**2*(beta2*X(1)+X(2));
  MAT(3,2)=-2*kf2*ima*pp*X(1);
  MAT(4,2)=(lambda2+M2*beta2**2)/Vf2**2*X(1)+2*mu2*kf2**2*X(1)+M2*beta2/Vf2**2*X(2);

  X=P2(:,2);
  MAT(1,3)=ks2*X(2)+ks2*X(1);
  MAT(2,3)=M2/Vs2**2*(beta2*X(1)+X(2));
  MAT(3,3)=-2*ks2*ima*pp*X(1);
  MAT(4,3)=(lambda2+M2*beta2**2)/Vs2**2*X(1)+2*mu2*ks2**2*X(1)+M2*beta2/Vs2**2*X(2);
  
  MAT(1,4)=ima*pp*(1.-rhof2/rhow2);
  MAT(2,4)=0;
  MAT(3,4)=(kpsi2**2+pp**2);
  MAT(4,4)=2*mu2*kpsi2*ima*pp;
  
  Coeff(1)=-1/(2*V1**2*rho1);
  Coeff(2)=-1/(2*k1*V1**2);
  Coeff(3)=0;
  Coeff(4)=-1/(2*k1*V1**2);

else
!!!! Sealed pore transmission conditions
   MAT(1,1)=-k1/rho1
   MAT(2,1)=0
   MAT(3,1)=0
   MAT(4,1)=1
   X=P2(:,1)
   MAT(1,2)=kf2*X(2)+kf2*X(1);
   MAT(2,2)=kf2*X(2)
   MAT(3,2)=-2*kf2*ima*pp*X(1);
   MAT(4,2)=(lambda2+M2*beta2**2)/Vf2**2*X(1)+2*mu2*kf2**2*X(1)+M2*beta2/Vf2**2*X(2);
   
  X=P2(:,2);
  MAT(1,3)=ks2*X(2)+ks2*X(1);
  MAT(2,3)=ks2*X(2)
  MAT(3,3)=-2*ks2*ima*pp*X(1);
  MAT(4,3)=(lambda2+M2*beta2**2)/Vs2**2*X(1)+2*mu2*ks2**2*X(1)+M2*beta2/Vs2**2*X(2);
  
  MAT(1,4)=ima*pp*(1.-rhof2/rhow2);
  MAT(2,4)=-ima*pp*rhof2/rhow2;
  MAT(3,4)=(kpsi2**2+pp**2);
  MAT(4,4)=2*mu2*kpsi2*ima*pp;
  
  Coeff(1)=-1/(2*V1**2*rho1);
  Coeff(2)=0
  Coeff(3)=0;
  Coeff(4)=-1/(2*k1*V1**2);
endif

  CALL ZGESV( 4, 1, MAT, 4, IPIV, Coeff, 4, INFO )
  calccoeff_ap=Coeff
end function calccoeff_ap
