c*********************************************************
c                  TF_param.inc                         **
c*********************************************************
      include 'phasespace.inc'
c
c external functions
c
        double precision pt,eta,phi,theta,rap,rho
	external pt,eta,phi,theta,rap,rho
c
c local
c
	double precision prov1,prov2,prov3,prov4,prov5
	double precision prov6,prov7,prov8,prov9,prov10
c
c	MET tag
c
	integer met_lhco,opt_lhco
	common/LHCO_met_tag/met_lhco,opt_lhco
c
c	generated
c
c
c     global
c
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)    
      common /to_diagram_kin/ momenta, mvir2
c
c parameters to init in transfert_card.dat
c	
        double precision  tf_lightjet_E_19
        double precision  tf_lightjet_E_20
        double precision  tf_lightjet_E_21
        double precision  tf_lightjet_E_22
        double precision  tf_lightjet_E_23
        double precision  tf_lightjet_E_24
        double precision  tf_lightjet_E_25
        double precision  tf_lightjet_E_26
        double precision  tf_lightjet_E_27
        double precision  tf_lightjet_E_28
        double precision  tf_lightjet_E_29
        double precision  tf_lightjet_E_30
        double precision  tf_lightjet_E_31
        double precision  tf_lightjet_E_32
        double precision  tf_lightjet_E_33
        double precision  tf_lightjet_E_34
        double precision  tf_lightjet_E_35
        double precision  tf_lightjet_E_36
        double precision  tf_bjet_E_1
        double precision  tf_bjet_E_2
        double precision  tf_bjet_E_3
        double precision  tf_bjet_E_4
        double precision  tf_bjet_E_5
        double precision  tf_bjet_E_6
        double precision  tf_bjet_E_7
        double precision  tf_bjet_E_8
        double precision  tf_bjet_E_9
        double precision  tf_bjet_E_10
        double precision  tf_bjet_E_11
        double precision  tf_bjet_E_12
        double precision  tf_bjet_E_13
        double precision  tf_bjet_E_14
        double precision  tf_bjet_E_15
        double precision  tf_bjet_E_16
        double precision  tf_bjet_E_17
        double precision  tf_bjet_E_18
       Common/to_TF_param/tf_lightjet_E_19,tf_lightjet_E_20,
     &tf_lightjet_E_21,tf_lightjet_E_22,tf_lightjet_E_23,
     &tf_lightjet_E_24,tf_lightjet_E_25,tf_lightjet_E_26,
     &tf_lightjet_E_27,tf_lightjet_E_28,tf_lightjet_E_29,
     &tf_lightjet_E_30,tf_lightjet_E_31,tf_lightjet_E_32,
     &tf_lightjet_E_33,tf_lightjet_E_34,tf_lightjet_E_35,
     &tf_lightjet_E_36,tf_bjet_E_1,tf_bjet_E_2,tf_bjet_E_3,
     &tf_bjet_E_4,tf_bjet_E_5,tf_bjet_E_6,tf_bjet_E_7,
     &tf_bjet_E_8,tf_bjet_E_9,tf_bjet_E_10,tf_bjet_E_11,
     &tf_bjet_E_12,tf_bjet_E_13,tf_bjet_E_14,tf_bjet_E_15,
     &tf_bjet_E_16,tf_bjet_E_17,tf_bjet_E_18