SUBROUTINE OUTPUT_ClCdCms
    USE Global
    IMPLICIT NONE
    
    INTEGER:: iBlock
    
    DO iBlock = 1, Max_Block
        CALL GLOBAL_SetPointers(iBlock)
        CALL OUTPUT_ClCdCm
    ENDDO

END SUBROUTINE OUTPUT_ClCdCms


SUBROUTINE OUTPUT_ClCdCm
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,M,N
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: NI2, NJ2, NK2

    INTEGER:: Ibgn,Iend,Jbgn,Jend,Kbgn,Kend
    INTEGER:: IJK,MinorMax

    INTEGER:: IJK_sur
    INTEGER::  KZ,ii
    REAL:: Timee
    REAL:: Clift,Cdrag,Cyaw,Cmx,Cmy,Cmz,Clp,Clv,Cdp,Cdv,Cz
    REAL:: Cfx,Cfy,Cfz

    REAL:: Sgn,Coff

    REAL:: aaa,bbb
    REAL:: P2,Cp,TT,RmiuW

    REAL:: cosa,sina,cosb,sinb
    REAL:: Xcs,Ycs,Zcs
    REAL:: Sx,Sy,Sz,St,Sl
    REAL:: Fnx_p,Fny_p,Fnz_p
    REAL:: dCl_pl,dCd_pl

    REAL:: Dn,Denw,Vcon
    REAL:: Sig,tau
    REAL:: dQ_x,dQ_y,dQ_z
    REAL:: Fnx_v, Fny_v, Fnz_v
    REAL:: dCl_vl, dCd_vl
    REAL:: dCx,dCy,dCz

    INTEGER:: Iw,Iwd,Iwall
    INTEGER:: Jw,Jwd,Jwall
    INTEGER:: Kw,Kwd,Kwall

	character(LEN=100):: FL_CLCDCM

    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1
    NI2=ThisBlock%NI2
    NJ2=ThisBlock%NJ2
    NK2=ThisBlock%NK2


	 DO N=1,ThisBlock%Num_wall
	     Ibgn=ThisBlock%IJK_BE(N,1)
	     Iend=ThisBlock%IJK_BE(N,2)
	     Jbgn=ThisBlock%IJK_BE(N,3)
	     Jend=ThisBlock%IJK_BE(N,4)
	     Kbgn=ThisBlock%IJK_BE(N,5)
	     Kend=ThisBlock%IJK_BE(N,6)
!
         IJK_sur=ThisBlock%IJK_sur(N)

	     if(IJK_sur == 1 ) IJK      = 1
	     if(IJK_sur == 1 ) MinorMax = 0
	     if(IJK_sur == 2 ) IJK      = 1
	     if(IJK_sur == 2 ) MinorMax = 1

	     if(IJK_sur == 3 ) IJK      = 2
	     if(IJK_sur == 3 ) MinorMax = 0
	     if(IJK_sur == 4 ) IJK      = 2
	     if(IJK_sur == 4 ) MinorMax = 1

	     if(IJK_sur == 5 ) IJK      = 3
	     if(IJK_sur == 5 ) MinorMax = 0
	     if(IJK_sur == 6 ) IJK      = 3
	     if(IJK_sur == 6 ) MinorMax = 1

!	     CALL cldm_blk001(Ibgn,Iend,Jbgn,Jend,Kbgn,Kend,IJK,MinorMax,
!     >          Clift,Cdrag,Cmx,Cmy,Cmz,Clp,Clv,Cdp,Cdv,Cz)


!!!!!!
!	  output subroutine for the Lift, drag and Pitching, Yawing and Roll moment
!!!!!!
        cosa  = cos(alfa)
        sina  = sin(alfa)
        cosb  = cos(beta_yaw)
        sinb  = sin(beta_yaw)
!!!!!!
!	  force in 3 direction
!!!!!!
	  Cfx  = 0.
	  Cfy  = 0.
	  Cfz  = 0.

	  Clp  = 0.
	  Clv  = 0.

	  Cdp  = 0.
	  Cdv  = 0.

	  Cmx  = 0.
	  Cmy  = 0.
	  Cmz  = 0.
!!!!!!
!	  Note: if (I==1  .or. J==1  .or. K==1  ) Sgn=-1.
!		    if (I==NI .or. J==NJ .or. K==NK ) Sgn=+1.
!!!!!!
	  if( MinorMax == 0)then	!! ±ß½çÎª1
	     Sgn  = -1.
	  else
	     Sgn  = +1.
	  endif
!!!!!!
	  St = 0.
!
	  DO K=Kbgn,Kend 
	  DO J=Jbgn,Jend 
	  DO I=Ibgn,Iend 
!!!!!!
!		 I direction
!!!!!!
	     if (IJK == 1) then
	        if(MinorMax == 0)then
	           Iw   = 1
	           Iwall= 1			!! pressure
	        else
	           Iw   = NI1
	           Iwall= NI		!! pressure
	        endif
	        Jw    = J
	        Jwall = J

	        Kw    = K
	        Kwall = K

           Xcs=(XX(Iw,J,K)+XX(Iw,J+1,K)+XX(Iw,J,K+1)+XX(Iw,J+1,K+1))/4.
           Ycs=(YY(Iw,J,K)+YY(Iw,J+1,K)+YY(Iw,J,K+1)+YY(Iw,J+1,K+1))/4.
           Zcs=(ZZ(Iw,J,K)+ZZ(Iw,J+1,K)+ZZ(Iw,J,K+1)+ZZ(Iw,J+1,K+1))/4.
	     endif
!!!!!!
!		 J direction
!!!!!!
	     if (IJK == 2) then
	        if(MinorMax == 0)then
	           Jw   = 1
	           Jwall= 1			!! pressure
	        else
	           Jw   = NJ1
	           Jwall= NJ		!! pressure
	        endif
	        Iw    = I
	        Iwall = I

	        Kw    = K
	        Kwall = K

           Xcs=(XX(I,Jw,K)+XX(I+1,Jw,K)+XX(I,Jw,K+1)+XX(I+1,Jw,K+1))/4.
           Ycs=(YY(I,Jw,K)+YY(I+1,Jw,K)+YY(I,Jw,K+1)+YY(I+1,Jw,K+1))/4.
           Zcs=(ZZ(I,Jw,K)+ZZ(I+1,Jw,K)+ZZ(I,Jw,K+1)+ZZ(I+1,Jw,K+1))/4.
	     endif
!!!!!!
!		 K direction
!!!!!!
	     if (IJK == 3) then
	        if(MinorMax == 0)then
	           Kw   = 1
	           Kwall= 1			!! pressure
	        else
	           Kw   = NK1
	           Kwall= NK		!! pressure
	        endif
	        Iw    = I
	        Iwall = I

	        Jw    = J
	        Jwall = J

           Xcs=(XX(I,J,Kw)+XX(I+1,J,Kw)+XX(I,J+1,Kw)+XX(I+1,J+1,Kw))/4.
           Ycs=(YY(I,J,Kw)+YY(I+1,J,Kw)+YY(I,J+1,Kw)+YY(I+1,J+1,Kw))/4.
           Zcs=(ZZ(I,J,Kw)+ZZ(I+1,J,Kw)+ZZ(I,J+1,Kw)+ZZ(I+1,J+1,Kw))/4.
	     endif
!!!!!!
!		 geometry
!!!!!!
		 Sx=SD(IJK,1,Iwall,Jwall,Kwall)/Grad(IJK,Iwall,Jwall,Kwall)
		 Sy=SD(IJK,2,Iwall,Jwall,Kwall)/Grad(IJK,Iwall,Jwall,Kwall)
		 Sz=SD(IJK,3,Iwall,Jwall,Kwall)/Grad(IJK,Iwall,Jwall,Kwall)
!		 St=SD(IJK,4,Iwall,Jwall,Kwall)/Grad(IJK,Iwall,Jwall,Kwall)
	     Sl=Grad(IJK,Iwall,Jwall,Kwall)
!!!!!!
!	     pressure
!!!!!!
	     P2     = PP(Iw,Jw,Kw)-PPF
           Fnx_p  =  P2 * Sx*Sl * Sgn
           Fny_p  =  P2 * Sy*Sl * Sgn
           Fnz_p  =  P2 * Sz*Sl * Sgn
!
!		 dCl_pl = Fny_p*cosa      - Fnx_p*sina
           dCd_pl = Fnx_p*cosa*cosb + Fny_p*sina*cosb + Fnz_p*sinb

           
!!!!!!
!	     Skin friction
!!!!!!
           TT     = T(Iw,Jw,Kw)
           RmiuW  = (1.+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5
!
		 Dn  = Dst(Iw,Jw,Kw)
		 Denw= V(1,Iw,Jw,Kw)
	     Vcon=(V(2,Iw,Jw,Kw)/Denw * Sx*Sl   &
     &          +V(3,Iw,Jw,Kw)/Denw * Sy*Sl &
     &          +V(4,Iw,Jw,Kw)/Denw * Sz*Sl)            ! * Sgn
!     &          +1.0000000000000000 * St*Sl)     ! * Sgn
!
	     Sig   = 1. !! sign(1., V(2,Iw,Jw,Kw))
	     tau   = Sig*RmiuW/Dn/Ref

		 dQ_x  = V(2,Iw,Jw,Kw)/Denw * Sl - Vcon * Sx !*Sl * sgn
		 dQ_y  = V(3,Iw,Jw,Kw)/Denw * Sl - Vcon * Sy !*Sl * sgn
		 dQ_z  = V(4,Iw,Jw,Kw)/Denw * Sl - Vcon * Sz !*Sl * sgn
!
           Fnx_v =  tau*dQ_x
           Fny_v =  tau*dQ_y
           Fnz_v =  tau*dQ_z
!
           dCd_vl= Fnx_v*cosa*cosb + Fny_v*sina*cosb + Fnz_v*sinb
           
           
!!!!!!
!		 integration for the force and moment
!!!!!!
	     dCx   =       Fnx_p  + Fnx_v
	     dCy   =       Fny_p  + Fny_v
	     dCz   =       Fnz_p  + Fnz_v
!
		 Cfx   = Cfx + Fnx_p  + Fnx_v
		 Cfy   = Cfy + Fny_p  + Fny_v
		 Cfz   = Cfz + Fnz_p  + Fnz_v
!
		 Cdp   = Cdp + dCd_pl					!!Pressure contribution for the drag coefficients
		 Cdv   = Cdv + dCd_vl
!
           Cmx   = Cmx + dCy*(Zcs-Zwc) - dCz*(Ycs-Ywc)
           Cmz   = Cmz + dCx*(Ycs-Ywc) - dCy*(Xcs-Xwc)
           Cmy   = Cmy + dCz*(Xcs-Xwc) - dCx*(Zcs-Zwc)
	  ENDDO
	  ENDDO
	  ENDDO
	  Cfx   = Cfx * 2.0_8/Area_ref
	  Cfy   = Cfy * 2.0_8/Area_ref
	  Cfz   = Cfz * 2.0_8/Area_ref
!
	  Cdp   = Cdp * 2.0_8/Area_ref
	  Cdv   = Cdv * 2.0_8/Area_ref
!
        Clift = Cfy*cosa      - Cfx*sina
        Cdrag =  Cfx*cosa*cosb + Cfy*sina*cosb +  Cfz*sinb
        Cyaw=   -Cfx*cosa*sinb - Cfy*sina*sinb +  Cfz*cosb
!
	  Cmx   = -Cmx *2.0_8/Area_ref/B_ref
	  Cmy   = -Cmy *2.0_8/Area_ref/B_ref
	  Cmz   =  Cmz *2.0_8/Area_ref/C_ref
     
KZ=MOD(KI,KQ)
      ThisBlock%CLDM(N,1,KZ)=Clift
      ThisBlock%CLDM(N,2,KZ)=Cdrag
      ThisBlock%CLDM(N,3,KZ)=Cyaw
      ThisBlock%CLDM(N,4,KZ)=Cmx
      ThisBlock%CLDM(N,5,KZ)=Cmy
      ThisBlock%CLDM(N,6,KZ)=Cmz
      ThisBlock%CLDM(N,7,KZ)=Cdp
      ThisBlock%CLDM(N,8,KZ)=Cdv

!!!!!! 

 
 IF(Kz  == 0 )THEN
         FL_CLCDCM = ThisBlock%FL_cldm (N)
!
 	     open(UNIT=18,FILE=trim(FilePathPrefix) // FL_CLCDCM,MODE='WRITE',ACCESS='APPEND', STATUS='UNKNOWN')  !,SHARE='DENYNONE')
!!!!!!
	 do ii=1,KQ-1
         timee=DT_LUSGS_Main*float(KI-KQ+ii)
	     WRITE(18,100)KI-KQ+ii,timee,ThisBlock%CLDM(N,1,ii),ThisBlock%CLDM(N,2,ii),ThisBlock%CLDM(N,3,ii),ThisBlock%CLDM(N,4,ii),ThisBlock%CLDM(N,5,ii),ThisBlock%CLDM(N,6,ii),ThisBlock%CLDM(N,7,ii),ThisBlock%CLDM(N,8,ii)
    enddo

         timee=DT_LUSGS_Main*float(KI)
	     WRITE(18,100)KI,timee,ThisBlock%CLDM(N,1,0),ThisBlock%CLDM(N,2,0),ThisBlock%CLDM(N,3,0),ThisBlock%CLDM(N,4,0),ThisBlock%CLDM(N,5,0),ThisBlock%CLDM(N,6,0),ThisBlock%CLDM(N,7,0),ThisBlock%CLDM(N,8,0)
	     close(18)
ENDIF
	  
      ENDDO
!!!!!
100	  format(1x,I6,9(2x,E16.8))

END SUBROUTINE OUTPUT_ClCdCm
