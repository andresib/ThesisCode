
C+-----------------------------------------------------------------------+
C|                  TRANSFER FUNCTION FOR MADWEIGHT                      |
C|                                                                       |
C|     Author: Pierre Artoisenet (UCL-CP3)                               |
C|             Olivier Mattelaer (UCL-CP3)                               |
C+-----------------------------------------------------------------------+
C|     This file is generated automaticly by MADWEIGHT-TF_BUILDER        |
C+-----------------------------------------------------------------------+     


C+-----------------------------------------------------------------------+
C|    Transfer function for tf_E_lightjet
C+-----------------------------------------------------------------------+
      subroutine tf_E_lightjet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'


        prov1=(tf_lightjet_E_19+tf_lightjet_E_20*dsqrt(p(0))+tf_lightjet_E_21*p(0))
        prov2=(tf_lightjet_E_22+tf_lightjet_E_23*dsqrt(p(0))+tf_lightjet_E_24*p(0))
        prov3=(tf_lightjet_E_25+tf_lightjet_E_26*dsqrt(p(0))+tf_lightjet_E_27*p(0))
        prov4=(tf_lightjet_E_28+tf_lightjet_E_29*dsqrt(p(0))+tf_lightjet_E_30*p(0))
        prov5=(tf_lightjet_E_31+tf_lightjet_E_32*dsqrt(p(0))+tf_lightjet_E_33*p(0))
        prov6=(tf_lightjet_E_34+tf_lightjet_E_35*dsqrt(p(0))+tf_lightjet_E_36*p(0))

        tf=(prov1*exp(-(p(0)-pexp(0)-prov2)**2/2d0/prov3**2))      !first gaussian
        tf=tf+(prov4*exp(-(p(0)-pexp(0)-prov5)**2/2d0/prov6**2))   !second gaussian
        tf=tf*((1d0/dsqrt(2d0*pi))/(prov1*prov3+prov4*prov6))      !normalisation



      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_E_lightjet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_E_lightjet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'


        prov3=(tf_lightjet_E_25+tf_lightjet_E_26*dsqrt(pexp(0))+tf_lightjet_E_27*pexp(0))
        prov6=(tf_lightjet_E_34+tf_lightjet_E_35*dsqrt(pexp(0))+tf_lightjet_E_36*pexp(0))

        width=max(prov3,prov6)



      width_E_lightjet= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_THETA_lightjet
C+-----------------------------------------------------------------------+
      subroutine tf_THETA_lightjet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_THETA_lightjet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_THETA_lightjet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'

        width=0d0


      width_THETA_lightjet= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_PHI_lightjet
C+-----------------------------------------------------------------------+
      subroutine tf_PHI_lightjet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_PHI_lightjet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_PHI_lightjet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'

        width=0d0


      width_PHI_lightjet= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_E_bjet
C+-----------------------------------------------------------------------+
      subroutine tf_E_bjet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'


        prov1=(tf_bjet_E_1+tf_bjet_E_2*dsqrt(p(0))+tf_bjet_E_3*p(0))
        prov2=(tf_bjet_E_4+tf_bjet_E_5*dsqrt(p(0))+tf_bjet_E_6*p(0))
        prov3=(tf_bjet_E_7+tf_bjet_E_8*dsqrt(p(0))+tf_bjet_E_9*p(0))
        prov4=(tf_bjet_E_10+tf_bjet_E_11*dsqrt(p(0))+tf_bjet_E_12*p(0))
        prov5=(tf_bjet_E_13+tf_bjet_E_14*dsqrt(p(0))+tf_bjet_E_15*p(0))
        prov6=(tf_bjet_E_16+tf_bjet_E_17*dsqrt(p(0))+tf_bjet_E_18*p(0))

        tf=(prov1*exp(-(p(0)-pexp(0)-prov2)**2/2d0/prov3**2))      !first gaussian
        tf=tf+(prov4*exp(-(p(0)-pexp(0)-prov5)**2/2d0/prov6**2))   !second gaussian
        tf=tf*((1d0/dsqrt(2d0*pi))/(prov1*prov3+prov4*prov6))      !normalisation 	



      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_E_bjet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_E_bjet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'


        prov3=(tf_bjet_E_4+tf_bjet_E_5*dsqrt(pexp(0))+tf_bjet_E_6*pexp(0))
        prov6=(tf_bjet_E_13+tf_bjet_E_14*dsqrt(pexp(0))+tf_bjet_E_15*pexp(0))

        width=max(prov3,prov6) 	



      width_E_bjet= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_THETA_bjet
C+-----------------------------------------------------------------------+
      subroutine tf_THETA_bjet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_THETA_bjet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_THETA_bjet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'

        width=0d0


      width_THETA_bjet= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_PHI_bjet
C+-----------------------------------------------------------------------+
      subroutine tf_PHI_bjet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_PHI_bjet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_PHI_bjet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'

        width=0d0


      width_PHI_bjet= width

      return
      end



