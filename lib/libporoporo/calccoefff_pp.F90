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

FUNCTION calccoefff_pp(pp)
  USE m_const
  USE m_phys
  IMPLICIT NONE
  COMPLEX*16,INTENT(in) :: pp
  INTEGER ::Info,IPIV(6)
  COMPLEX*16 :: kf1,ks1,kpsi1,kf2,ks2,kpsi2
  COMPLEX*16 :: Mat(6,6),coeff(6),calccoefff_pp(6),X(2)
  IF(REAL(pp).NE.0.d0) THEN
     kf1=SQRT(pp**2+1/Vf1**2)
     ks1=SQRT(pp**2+1/Vs1**2)
     kpsi1=SQRT(pp**2+1/Vpsi1**2)
     kf2=SQRT(pp**2+1/Vf2**2)
     ks2=SQRT(pp**2+1/Vs2**2)
     kpsi2=SQRT(pp**2+1/Vpsi2**2)
  ELSE
     IF (AIMAG(pp)**2>1/Vf1**2) THEN
        kf1=ima*SQRT(-pp**2-1/Vf1**2)
     ELSE
        kf1=SQRT(pp**2+1/Vf1**2)
     END IF
     IF (AIMAG(pp)**2>1/Vpsi1**2) THEN
        kpsi1=ima*SQRT(-pp**2-1/Vpsi1**2)
     ELSE
        kpsi1=SQRT(pp**2+1/Vpsi1**2)
     END IF
     IF (AIMAG(pp)**2>1/Vs1**2) THEN
        ks1=ima*SQRT(-pp**2-1/Vs1**2)
     ELSE
        ks1=SQRT(pp**2+1/Vs1**2)
     END IF
     IF (AIMAG(pp)**2>1/Vf2**2) THEN
        kf2=ima*SQRT(-pp**2-1/Vf2**2)
     ELSE
        kf2=SQRT(pp**2+1/Vf2**2)
     END IF
     IF (AIMAG(pp)**2>1/Vs2**2) THEN
        ks2=ima*SQRT(-pp**2-1/Vs2**2)
     ELSE
        ks2=SQRT(pp**2+1/Vs2**2)
     END IF
     IF (AIMAG(pp)**2>1/Vpsi2**2) THEN
        kpsi2=ima*SQRT(-pp**2-1/Vpsi2**2)
     ELSE
        kpsi2=SQRT(pp**2+1/Vpsi2**2)
     END IF
  END IF
  X=P1(:,1)
  MAT(1,1)=-ima*pp*X(1)
  MAT(2,1)=-kf1*X(1)
  MAT(3,1)=-kf1*X(2)
  MAT(4,1)=M1/Vf1**2*(beta1*X(1)+X(2))
  MAT(5,1)=2*mu1*kf1*ima*pp*X(1)
  MAT(6,1)=(lambda1+M1*beta1**2)/Vf1**2*X(1)+2*mu1*kf1**2*X(1)+M1*beta1/Vf1**2*X(2);

  X=P1(:,2)
  MAT(1,2)=-ima*pp*X(1)
  MAT(2,2)=-ks1*X(1)
  MAT(3,2)=-ks1*X(2)
  MAT(4,2)=M1/Vs1**2*(beta1*X(1)+X(2))
  MAT(5,2)=2*mu1*ks1*ima*pp*X(1)
  MAT(6,2)=(lambda1+M1*beta1**2)/Vs1**2*X(1)+2*mu1*ks1**2*X(1)+M1*beta1/Vs1**2*X(2)

  MAT(1,3)=-kpsi1
  MAT(2,3)=ima*pp
  MAT(3,3)=-ima*pp*rhof1/rhow1
  MAT(4,3)=0
  MAT(5,3)=mu1*(kpsi1**2+pp**2)
  MAT(6,3)=-2*mu1*kpsi1*ima*pp

  X=P2(:,1)
  MAT(1,4)=ima*pp*X(1)
  MAT(2,4)=-kf2*X(1)
  MAT(3,4)=-kf2*X(2)
  MAT(4,4)=-M2/Vf2**2*(beta2*X(1)+X(2))
  MAT(5,4)=2*mu2*kf2*ima*pp*X(1)
  MAT(6,4)=-(lambda2+M2*beta2**2)/Vf2**2*X(1)-2*mu2*kf2**2*X(1)-M2*beta2/Vf2**2*X(2)
  X=P2(:,2)
  MAT(1,5)=ima*pp*X(1)
  MAT(2,5)=-ks2*X(1)
  MAT(3,5)=-ks2*X(2)
  MAT(4,5)=-M2/Vs2**2*(beta2*X(1)+X(2))
  MAT(5,5)=2*mu2*ks2*ima*pp*X(1)
  MAT(6,5)=-(lambda2+M2*beta2**2)/Vs2**2*X(1)-2*mu2*ks2**2*X(1)-M2*beta2/Vs2**2*X(2)

  MAT(1,6)=-kpsi2
  MAT(2,6)=-ima*pp
  MAT(3,6)=ima*pp*rhof2/rhow2
  MAT(4,6)=0
  MAT(5,6)=-mu2*(kpsi2**2+pp**2)
  MAT(6,6)=-2*mu2*kpsi2*ima*pp

  X=P1(:,1)/(2*kf1)
  Coeff(1)=ima*pp*X(1)
  Coeff(2)=-kf1*X(1)
  Coeff(3)=-kf1*X(2)
  Coeff(4)=-M1/Vf1**2*(beta1*X(1)+X(2))
  Coeff(5)=2*mu1*ima*pp*kf1*X(1)
  Coeff(6)=-(lambda1+M1*beta1**2)/Vf1**2*X(1)-2*mu1*kf1**2*X(1)-M1*beta1/Vf1**2*X(2)
  Coeff=Coeff/Vf1**2
  CALL ZGESV( 6, 1, MAT, 6, IPIV, Coeff, 6, INFO )
  calccoefff_pp=Coeff
END FUNCTION calccoefff_pp
