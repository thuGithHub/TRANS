SUBROUTINE TIME_SpectralRadius
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: KKtm
    INTEGER:: II,JJ
    
    REAL:: CFL_max, CFL_min
    REAL:: dT_max, dT_min
    REAL:: Den,Vxi,Vyi,Vzi
    REAL:: AAc,UU,VV,WW
    REAL:: TT
    REAL:: YY1,YY2,YY3,YY4,YY5,YY6
    REAL:: RmiuL,RmiuTL,uvmu,uvmupr,pruvmu
    REAL:: ARX1,ARY1,ARZ1
    REAL:: ARX2,ARY2,ARZ2
    REAL:: ARX3,ARY3,ARZ3
    REAL:: ST1sq, ST2sq, ST3sq
    REAL:: UVW2,ubetac2,ALL,ubetac2u,ubetac2v,ubetac2w,RmiuAll
    REAL:: alfu,AACu
    REAL:: alfv,AACv
    REAL:: alfw,AACw
            
    REAL:: ST1,ST11,ST2,ST22,ST3,ST33
    REAL:: Dtm_local, CFL_local
    real::vibc,vjbc,vkbc

    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1

!!!!!!
!	  spectral radius
!!!!!!
	  CFL_max=0.
	  CFL_min=1.E20
	   dT_max=0.
	   dT_min=1.E20
!!!!!!
!        DO K=1,NK1
        DO K=1,NK1
        DO J=1,NJ1
        DO I=1,NI1
!
           Den=V(1,I,J,K)
           Vxi=V(2,I,J,K)/Den
           Vyi=V(3,I,J,K)/Den
           Vzi=V(4,I,J,K)/Den
           !AAc =SQRT(gam(I,J,K)*PP(I,J,K)/Den)
           AAc =SQRT(1.4*PP(I,J,K)/Den)  !by ydd

		 TT = T(I,J,K)
         RmiuL = (1.0+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5
         
		 RmiuTL = Rmiu(I,J,K)
		 uvmu = (RmiuL + RmiuTL + 1.0e-20)*4.0/3.0
	     !uvmupr = (RmiuL/0.72 + RmiuTL/0.9 + 1.0e-20)*gam(I,J,K)
	     uvmupr = (RmiuL/0.72 + RmiuTL/0.9 + 1.0e-20)*1.4   !by ydd
		 pruvmu = max(uvmu, uvmupr) / Den

		 ARX1 = 0.5*(SD(1,1,I,J,K)+SD(1,1,I+1,J,K))
		 ARY1 = 0.5*(SD(1,2,I,J,K)+SD(1,2,I+1,J,K))
		 ARZ1 = 0.5*(SD(1,3,I,J,K)+SD(1,3,I+1,J,K))
		 ST1sq  = ARX1*ARX1+ARY1*ARY1+ARZ1*ARZ1
!added by ydd
        vibc=0.5*(vibn(i,j,k)+vibn(i+1,j,k))

		 ARX2 = 0.5*(SD(2,1,I,J,K)+SD(2,1,I,J+1,K))
		 ARY2 = 0.5*(SD(2,2,I,J,K)+SD(2,2,I,J+1,K))
		 ARZ2 = 0.5*(SD(2,3,I,J,K)+SD(2,3,I,J+1,K))
		 ST2sq  = ARX2*ARX2+ARY2*ARY2+ARZ2*ARZ2
        vjbc=0.5*(vjbn(i,j,k)+vjbn(i,j+1,k))

		 ARX3 = 0.5*(SD(3,1,I,J,K)+SD(3,1,I,J,K+1))
		 ARY3 = 0.5*(SD(3,2,I,J,K)+SD(3,2,I,J,K+1))
		 ARZ3 = 0.5*(SD(3,3,I,J,K)+SD(3,3,I,J,K+1))
		 ST3sq  = ARX3*ARX3+ARY3*ARY3+ARZ3*ARZ3
        vkbc=0.5*(vkbn(i,j,k)+vkbn(i,j,k+1))

		 UU = Vxi*ARX1+Vyi*ARY1+Vzi*ARZ1+vibc
		 VV = Vxi*ARX2+Vyi*ARY2+Vzi*ARZ2+vjbc
		 WW = Vxi*ARX3+Vyi*ARY3+Vzi*ARZ3+vkbc


         IF(IF_PRECONDITION) THEN !IF_PRECONDITION) THEN !low Ma

		     UVW2 = Vxi*Vxi+Vyi*Vyi+Vzi*Vzi
	         ubetac2 = 0.0_8
		     ALL = min(Alagm(1,I,J,K),Alagm(2,I,J,K),Alagm(3,I,J,K))
		     CALL getPrecondPara(UVW2,AAc*AAc,Den, 0.0, 0.0,1.0, ubetac2)

		     ubetac2u = ubetac2
		     ubetac2v = ubetac2
		     ubetac2w = ubetac2

		     RmiuAll = RmiuL+ RmiuTL

		     ubetac2u = max(ubetac2, (RmiuAll/Den/Ref/Alagm(1,I,J,K)/AAc)**2.0)
		     alfu = (1.0-ubetac2u) / 2.0
		     AAcu=sqrt(alfu*alfu*UVW2+ubetac2u*AAc*AAc)

		     ubetac2v = max(ubetac2, (RmiuAll/Den/Ref/Alagm(2,I,J,K)/AAc)**2.0)
		     alfv = (1.0-ubetac2v) / 2.0
		     AAcv=sqrt(alfv*alfv*UVW2+ubetac2v*AAc*AAc)

		     ubetac2w = max(ubetac2, (RmiuAll/Den/Ref/Alagm(3,I,J,K)/AAc)**2.0)
		     alfw = (1.0-ubetac2w) / 2.0
		     AAcw=sqrt(alfw*alfw*UVW2+ubetac2w*AAc*AAc)

		ELSE
			alfu =0.0
			alfv =0.0
			alfw =0.0

			AAcu=AAc
			AAcv=AAc
			AAcw=AAc
		ENDIF




           Rds(1,I,J,K)=abs(UU)*(1.0_8-alfu) + sqrt(ST1sq)*AAcu +   &
     &	 	  2.0_8*ST1sq/Vol(I,J,K)*pruvmu / Ref
           Rds(2,I,J,K)=abs(VV)*(1.0_8-alfv) + sqrt(ST2sq)*AAcv +   &
     &	 	  2.0_8*ST2sq/Vol(I,J,K)*pruvmu / Ref
           Rds(3,I,J,K)=abs(WW)*(1.0_8-alfw) + sqrt(ST3sq)*AAcw +   &
     &	 	  2.0_8*ST3sq/Vol(I,J,K)*pruvmu / Ref


	     Dtm_local=Vol(I,J,K)/(Rds(1,I,J,K)+Rds(2,I,J,K)+Rds(3,I,J,K))


!!!!!!
!		local or global time step
!!!!!!
        

		
         if(IF_local_Main == 1 )then
		    Dtm(I,J,K) = CFL_Main * Dtm_local
	        
			if(Dtm(I,J,K) > dT_max) dT_max=Dtm(I,J,K)
			if(Dtm(I,J,K) < dT_min) dT_min=Dtm(I,J,K)
	     else
		    Dtm(I,J,K) = DT_LUSGS_Main
	        CFL_local  = DT_LUSGS_Main / Dtm_local

			if(CFL_local > CFL_max)CFL_max=CFL_local
			if(CFL_local < CFL_min)CFL_min=CFL_local
	     endif
!!!!!!
	     CFL_sub=CFL_sub_start +(CFL_sub_end- CFL_sub_start)*(KI-KI_ini)/CFL_sub_step
         CFL_sub=min(CFL_sub,CFL_sub_end)
         
         if (Kind_dual /= 0 ) then
		    if(IF_local_Sub  == 1 )then
		       
               Dtm(I,J,K) = CFL_Sub * Dtm_local
	        
			   if(Dtm(I,J,K) > dT_max) dT_max=Dtm(I,J,K)
			   if(Dtm(I,J,K) < dT_min) dT_min=Dtm(I,J,K)
	        else
		       Dtm(I,J,K) = DT_LUSGS_Sub 
	           CFL_local  = DT_LUSGS_Sub / Dtm_local

			   if(CFL_local > CFL_max)CFL_max=CFL_local
			   if(CFL_local < CFL_min)CFL_min=CFL_local
	        endif
	     endif
!!!!!!
	  ENDDO
	  ENDDO
	  ENDDO
!!!!!!
        KKtm=MOD(KI,KO)
!        IF(KKtm == 0) THEN
!	     if(If_local_main == 0 .or. If_local_sub == 0 )then
!             write(*,*)" CFL_max= ",CFL_max," CFL_min= ",CFL_min
!	     else
!        IF(Kchld ==1)THEN
!             write(*,*)"  Dt_max= ", dT_max,"  dT_min= ", dT_min, " CFL_local=",CFL_Sub
!          else
!             write(*,*)"  CFL_max= ", CFL_max,"  CFL_min= ", CFL_min, " DT_LUSGS =", DT_LUSGS_Main
!	     endif
!	  endif
!!!!!!
!        DO K=1,NK1      !spectral radius for the first ghost cells, by ydd
!        DO J=1,NJ1
!        DO I=0,NI,NI
!           Den=V(1,I,J,K)
!           Vxi=V(2,I,J,K)/Den
!           Vyi=V(3,I,J,K)/Den
!           Vzi=V(4,I,J,K)/Den
!           AAc =SQRT(1.4*PP(I,J,K)/Den)  !by ydd
!            TT = T(I,J,K)
!            RmiuL = (1.0+Csthlnd)/(TT+Csthlnd)*TT**1.5
!            RmiuTL = Rmiu(I,J,K)
!            uvmu = (RmiuL + RmiuTL + 1.0e-20)*4.0/3.0
            !uvmupr = (RmiuL/0.72 + RmiuTL/0.9 + 1.0e-20)*gam(I,J,K)
!            uvmupr = (RmiuL/0.72 + RmiuTL/0.9 + 1.0e-20)*1.4   !by ydd
!            pruvmu = max(uvmu, uvmupr) / Den

 !           ARX1 = 0.5*(SD(1,1,I,J,K)+SD(1,1,I+1,J,K))
!            ARY1 = 0.5*(SD(1,2,I,J,K)+SD(1,2,I+1,J,K))
!            ARZ1 = 0.5*(SD(1,3,I,J,K)+SD(1,3,I+1,J,K))
!            ST1sq  = ARX1*ARX1+ARY1*ARY1+ARZ1*ARZ1

!            ARX2 = 0.5*(SD(2,1,I,J,K)+SD(2,1,I,J+1,K))
!            ARY2 = 0.5*(SD(2,2,I,J,K)+SD(2,2,I,J+1,K))
!            ARZ2 = 0.5*(SD(2,3,I,J,K)+SD(2,3,I,J+1,K))
!            ST2sq  = ARX2*ARX2+ARY2*ARY2+ARZ2*ARZ2
!
!            ARX3 = 0.5*(SD(3,1,I,J,K)+SD(3,1,I,J,K+1))
!            ARY3 = 0.5*(SD(3,2,I,J,K)+SD(3,2,I,J,K+1))
!            ARZ3 = 0.5*(SD(3,3,I,J,K)+SD(3,3,I,J,K+1))
!            ST3sq  = ARX3*ARX3+ARY3*ARY3+ARZ3*ARZ3
!            vkbc=0.5*(vkbn(i,j,k)+vkbn(i,j,k+1))
!
!            UU = Vxi*ARX1+Vyi*ARY1+Vzi*ARZ1+vibc
!            VV = Vxi*ARX2+Vyi*ARY2+Vzi*ARZ2+vjbc
!            WW = Vxi*ARX3+Vyi*ARY3+Vzi*ARZ3+vkbc
!
!            IF(IF_PRECONDITION) THEN !IF_PRECONDITION) THEN !low Ma
!                UVW2 = Vxi*Vxi+Vyi*Vyi+Vzi*Vzi
!                ubetac2 = 0.0_8
!                ALL = min(Alagm(1,I,J,K),Alagm(2,I,J,K),Alagm(3,I,J,K))
!                CALL getPrecondPara(UVW2,AAc*AAc,Den, 0.0, 0.0,1.0, ubetac2!)

!                ubetac2u = ubetac2
!                ubetac2v = ubetac2
!                ubetac2w = ubetac2
!
!                RmiuAll = RmiuL+ RmiuTL
!
!                ubetac2u = max(ubetac2, (RmiuAll/Den/Ref/Alagm(1,I,J,K)/AAc)**2.0)
!                alfu = (1.0-ubetac2u) / 2.0
!                AAcu=sqrt(alfu*alfu*UVW2+ubetac2u*AAc*AAc)
!
!                ubetac2v = max(ubetac2, (RmiuAll/Den/Ref/Alagm(2,I,J,K)/AAc)**2.0)
!                alfv = (1.0-ubetac2v) / 2.0
 !               AAcv=sqrt(alfv*alfv*UVW2+ubetac2v*AAc*AAc)

!                ubetac2w = max(ubetac2, (RmiuAll/Den/Ref/Alagm(3,I,J,K)/AAc)**2.0)
!                alfw = (1.0-ubetac2w) / 2.0
!                AAcw=sqrt(alfw*alfw*UVW2+ubetac2w*AAc*AAc)
!
!            ELSE
!                alfu =0.0
!                alfv =0.0
!                alfw =0.0
!
!                AAcu=AAc
!                AAcv=AAc
!                AAcw=AAc
!            ENDIF
!
!           Rds(1,I,J,K)=abs(UU)*(1.0_8-alfu) + sqrt(ST1sq)*AAcu +   &
!     &        2.0_8*ST1sq/Vol(I,J,K)*pruvmu / Ref
!           Rds(2,I,J,K)=abs(VV)*(1.0_8-alfv) + sqrt(ST2sq)*AAcv +   &
!     &        2.0_8*ST2sq/Vol(I,J,K)*pruvmu / Ref
!           Rds(3,I,J,K)=abs(WW)*(1.0_8-alfw) + sqrt(ST3sq)*AAcw +   &
!     &        2.0_8*ST3sq/Vol(I,J,K)*pruvmu / Ref
!        enddo
!        enddo
!        enddo
!
END SUBROUTINE TIME_SpectralRadius
