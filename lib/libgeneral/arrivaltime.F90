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

FUNCTION arrivaltime(h,x,y,V1,V2)
IMPLICIT NONE
REAL*8,INTENT(in) ::h,x,y,V1,V2
REAL*8 :: arrivaltime,alpha,Rwork(8),zz(4)
COMPLEX*16 :: Mat(4,4),Eig(4),Work(100)
INTEGER ::Info,k,ii
Mat=0.D0
IF (x.NE.0) THEN
   IF (V1.NE.V2) THEN
      alpha=(1/(V1**2)-1/(V2**2))
      Mat(1,1)=-(2*x/V2**2-2/V1**2*x)/alpha
      Mat(1,2)=-(-x**2/V2**2-1/V2**2*h**2+1/V1**2*(x**2+y**2))/alpha
      Mat(1,3)=-(2*x/V2**2*h**2)/alpha
      Mat(1,4)=-(-x**2/V2**2*h**2)/alpha
      Mat(2,1)=1.D0
      Mat(3,2)=1.D0
      Mat(4,3)=1.D0
      CALL ZGEEV('N', 'N', 4, Mat, 4, Eig, 1, 1,1,1,&
           &WORK, 100, RWORK, INFO )
      ii=0
      DO k=1,4
         IF (AIMAG(eig(k)).EQ.0) THEN
            IF (x*REAL(eig(k))>0) THEN
               ii=ii+1
               zz(ii)=REAL(eig(k))
            END IF
         ELSEIF (REAL(eig(k)).NE.0) THEN
            IF (ABS(AIMAG(eig(k))/ABS(eig(k))).LT.1e-6) THEN
               IF (x*REAL(eig(k))>0) THEN
                  ii=ii+1
                  zz(ii)=REAL(eig(k))
               END IF
            END IF
         END IF
      END DO
   ELSE
      ii=1
      zz=(x*h)/(-y+h)
   END IF
ELSE
   ii=1
   zz=0.D0
END IF
arrivaltime=MINVAL(SQRT(h**2+zz(1:ii)**2)/V1+SQRT((x-zz(1:ii))**2+y**2)/V2)
END FUNCTION arrivaltime
