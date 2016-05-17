SUBROUTINE TIME_RK
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,M
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: II,JJ
    
    REAL:: DRS,DUa,DR,DRSa
    INTEGER:: Ima,Jma,Kma
    
    
    
!    tol=1.0E-4
    
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1
    
    DRSa=0.0
    DUa=0.0
    DO K=1,NK1
           DO J=1,NJ1
              DO I=1,NI1
                 
                 DR=ABS(D(1,I,J,K)/V(1,I,J,K))
                 
                 IF(ISNAN(D(1,I,J,K)))then
                     write(*,*) I,J,K
                     write(*,*) KI, D(1,I,J,K)
                     continue
                 endif
                 
                 IF(DUa<=DR)THEN
                    DUa=DR
                    Ima=I
                    Jma=J
                    Kma=K
	           ENDIF
                 DRSa=DRSa+DR*DR
	        ENDDO 
	     ENDDO
        ENDDO  

        DRSa=SQRT(DRSa/NI1/NJ1/NK1)

        if(ISNAN(DRSa)) then
        write(*,*)"Residual of blk",ThisBlock%id_present_blk,"is NAN"
        continue
        stop
        endif
    
    do K=1,NK1
	do J=1,NJ1
	do I=1,NI1
        
        IF(Itr_Chld == 1)then
            
            do L=1,5
                V(L,I,J,K)=V(L,I,J,K)+ Dtm(I,J,K)*D(L,I,J,K)/Vol(I,J,K)
            enddo
        
        ENDIF
        
        IF(Itr_Chld == 2)then
            do L=1,5
                V(L,I,J,K)=0.75*Wt2(L,I,J,K)+0.25*V(L,I,J,K)+ 0.25*Dtm(I,J,K)*D(L,I,J,K)/Vol(I,J,K)
            enddo
        ENDIF
        
        IF(Itr_Chld == 3)then
            do L=1,5
                V(L,I,J,K)=1.0/3.0*Wt2(L,I,J,K)+2.0/3.0*V(L,I,J,K)+ 2.0/3.0*Dtm(I,J,K)*D(L,I,J,K)/Vol(I,J,K)
            enddo
        ENDIF
        
    ENDDO
    ENDDO
    ENDDO
    
  
    
    ThisBlock%DRSa=DRSa
    ThisBlock%DUa=DUa
    ThisBlock%Ima=Ima
    ThisBlock%Jma=Jma
    ThisBlock%Kma=Kma
    
    
    END SUBROUTINE
