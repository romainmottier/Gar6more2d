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

function calccoefff_wall_pp(pp)
  Use m_const
  Use m_phys
implicit none
complex*16,intent(in) :: pp
integer ::Info,IPIV(3)
complex*16 :: kf1,ks1,kpsi1
complex*16 :: Mat(3,3),coeff(3),calccoefff_wall_pp(3),X(2)
if(real(pp).ne.0.d0) then
   kf1=sqrt(pp**2+1/Vf1**2)
   ks1=sqrt(pp**2+1/Vs1**2)
   kpsi1=sqrt(pp**2+1/Vpsi1**2)
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
 end if
 X=P1(:,1)
 MAT(1,1)=-ima*pp*X(1)
 MAT(2,1)=-kf1*X(1)
 MAT(3,1)=-kf1*X(2)

  X=P1(:,2)
  MAT(1,2)=-ima*pp*X(1)
  MAT(2,2)=-ks1*X(1)
  MAT(3,2)=-ks1*X(2)
  
  MAT(1,3)=-kpsi1
  MAT(2,3)=ima*pp
  MAT(3,3)=-ima*pp*rhof1/rhow1
  
  X=P1(:,1)/(2*kf1)
  Coeff(1)=ima*pp*X(1)
  Coeff(2)=-kf1*X(1)
  Coeff(3)=-kf1*X(2)
  Coeff=Coeff/Vf1**2
  CALL ZGESV( 3, 1, MAT, 3, IPIV, Coeff, 3, INFO )
  calccoefff_wall_pp=Coeff
end function calccoefff_wall_pp
