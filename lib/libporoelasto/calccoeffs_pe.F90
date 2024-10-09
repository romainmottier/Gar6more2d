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

function calccoeffs_pe(pp)
  use m_const
  use m_phys
  implicit none
  complex*16,intent(in) :: pp
  integer ::Info,IPIV(5)
  complex*16 :: kf1,ks1,kpsi1,kp2,kS2,X(2)
  complex*16 :: Mat(5,5),coeff(5),calccoeffs_pe(5)
  if(real(pp).ne.0.d0) then
     kf1=sqrt(pp**2+1/Vf1**2)
     ks1=sqrt(pp**2+1/Vs1**2)
     kpsi1=sqrt(pp**2+1/Vpsi1**2)
     kp2=sqrt(pp**2+1/Vp2**2)
     ks2=sqrt(pp**2+1/Vs2**2)
  else
     if (aimag(pp)**2>1/Vf1**2) then
        kf1=ima*sqrt(-pp**2-1/Vf1**2)
     else
        kf1=sqrt(pp**2+1/Vf1**2)
     end if
     if (aimag(pp)**2>1/Vpsi1**2) then
        kpsi1=ima*sqrt(-pp**2-1/Vpsi1**2)
     else
        kpsi1=sqrt(pp**2+1/Vpsi1**2)
     end if
     if (aimag(pp)**2>1/Vs1**2) then
        ks1=ima*sqrt(-pp**2-1/Vs1**2)
     else
        ks1=sqrt(pp**2+1/Vs1**2)
     end if
     if (aimag(pp)**2>1/Vf2**2) then
        kp2=ima*sqrt(-pp**2-1/Vp2**2)
     else
        kp2=sqrt(pp**2+1/Vp2**2)
     end if
     if (aimag(pp)**2>1/Vs2**2) then
        ks2=ima*sqrt(-pp**2-1/Vs2**2)
     else
        ks2=sqrt(pp**2+1/Vs2**2)
     end if
  end if
  X=P1(:,1)
  MAT(1,1)=-ima*pp*X(1)
  MAT(2,1)=-kf1*X(1)
  MAT(3,1)=kf1*X(2)!-M1/Vf1**2*(beta1*X(1)+X(2))!-kf1*X(2)
  MAT(4,1)=2*mu1*kf1*ima*pp*X(1)
  MAT(5,1)=((lambda1+M1*beta1**2)/Vf1**2*X(1)&
       &+2*mu1*kf1**2*X(1)+M1*beta1/Vf1**2*X(2));

  X=P1(:,2)
  MAT(1,2)=-ima*pp*X(1)
  MAT(2,2)=-ks1*X(1)
  MAT(3,2)=ks1*X(2)!-M1/Vs1**2*(beta1*X(1)+X(2))!-ks1*X(2)
  MAT(4,2)=2*mu1*ks1*ima*pp*X(1)
  MAT(5,2)=((lambda1+M1*beta1**2)/Vs1**2*X(1)+2*mu1*ks1**2*X(1)+M1*beta1/Vs1**2*X(2))

  MAT(1,3)=-kpsi1
  MAT(2,3)=ima*pp
  MAT(3,3)=ima*pp*rhof1/rhow1!0!-ima*pp*rhof1/rhow1
  MAT(4,3)=mu1*(kpsi1**2+pp**2)
  MAT(5,3)=-2*mu1*kpsi1*ima*pp

  MAT(1,4)=ima*pp
  MAT(2,4)=-kp2
  MAT(3,4)=0!-lambda2/Vp2**2-2*mu2*kp2**2
  MAT(4,4)=2*mu2*kp2*ima*pp
  MAT(5,4)=-lambda2/Vp2**2-2*mu2*kp2**2

  MAT(1,5)=-ks2
  MAT(2,5)=-ima*pp
  MAT(3,5)=0!-2*mu2*ks2*ima*pp
  MAT(4,5)=-mu2*(ks2**2+pp**2)
  MAT(5,5)=-2*mu2*ks2*ima*pp


  X=P1(:,2)/(2*ks1)
  Coeff(1)=ima*pp*X(1)
  Coeff(2)=-ks1*X(1)
  Coeff(3)=ks1*X(2)!-M1/Vs1**2*(beta1*X(1)+X(2))!-ks1*X(2)
  Coeff(4)=2*mu1*ima*pp*ks1*X(1)
  Coeff(5)=-((lambda1+M1*beta1**2)/Vs1**2*X(1)+2*mu1*ks1**2*X(1)+M1*beta1/Vs1**2*X(2))
  Coeff=Coeff/Vs1**2
  call ZGESV( 5, 1, MAT, 5, IPIV, Coeff, 5, INFO )
  calccoeffs_pe=Coeff
end function calccoeffs_pe
