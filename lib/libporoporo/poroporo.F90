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

SUBROUTINE poroporo
  Use m_phys  
  Use m_num
  Use m_source
  Use m_sismo
  implicit none 
        write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(6,*) '!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!'
        write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(6,*) 'In the poroelastics case, &
             &the code does not really compute the displacement, &
             &but its integral (for some reasons related to &
             &the Cagniard-de Hoop method, see the documentation). &
             &Therefore, don''t forget to replace the source by the derivative &
             &of the source function you are using.'
        write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!!! Computation of the incident wave
  allocate(Ux(1:nx,1:nt))
  allocate(Uy(1:nx,1:nt))
  Ux=0.d0
  Uy=0.d0
  if (type_medium.eq.1) then
!!!Infinite Medium
        CALL sub_incid_pp_f
        CALL sub_incid_pp_s
  elseif (type_medium.eq.2) then
!!!Free Surface
     if (y.ge.0d0) then
        CALL sub_incid_pp_f
        CALL sub_incid_pp_s
        CALL sub_reflexff_free_pp
        CALL sub_reflexfpsi_free_pp
        CALL sub_reflexfs_free_pp
        CALL sub_reflexsf_free_pp
        CALL sub_reflexspsi_free_pp
        CALL sub_reflexss_free_pp
     else
        write(6,*) 'Error, y can not be smaller than 0 in this configuration'
        stop
     end if
  elseif (type_medium.eq.3) then
!!!Wall Boundary
     if (y.ge.0d0) then
        CALL sub_incid_pp_f
        CALL sub_incid_pp_s
        CALL sub_reflexff_wall_pp
        CALL sub_reflexfpsi_wall_pp
        CALL sub_reflexfs_wall_pp
        CALL sub_reflexsf_wall_pp
        CALL sub_reflexspsi_wall_pp
        CALL sub_reflexss_wall_pp
     else
        write(6,*) 'Error, y can not be smaller than 0 in this configuration'
        stop
     end if
  elseif (type_medium.eq.4) then
     if (y.ge.0d0) then
        CALL sub_incid_pp_f
        CALL sub_incid_pp_s
        CALL sub_reflexff_pp
        CALL sub_reflexfpsi_pp
        CALL sub_reflexfs_pp
        CALL sub_reflexsf_pp
        CALL sub_reflexspsi_pp
        CALL sub_reflexss_pp
     else
        CALL sub_transmitff_pp
        CALL sub_transmitfpsi_pp
        CALL sub_transmitfs_pp
        CALL sub_transmitsf_pp
        CALL sub_transmitspsi_pp
        CALL sub_transmitss_pp
     end if
  end if
end SUBROUTINE poroporo
