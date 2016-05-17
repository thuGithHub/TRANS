SUBROUTINE COMMON_Define_BCinput_Parameter
    USE Global
    IMPLICIT NONE

    !	i_suf = 1 - 6: imin,imax,jmin,jmax,kmin,kmax
    !   i_Kind with i_Kindsub: Kind of BC
    !   i_Blk_Ex, i_Blkbd_Ex: p2p or quasip2p, correspondant domain BC
    !	  i_lcross 
    !	  i_lrev_a 
    !	  i_lrev_b 
                         
    I_Suf=1             !IJKBC(:,I_Suf),������ڶ��е�һ�д����߽�6��λ�����ͣ�imin,imax,jmin,jmax,kmin,kmax        

    I_IBgn=2            
    I_Iend=3            
    I_JBgn=4           
    I_Jend=5           
    I_KBgn=6            
    I_Kend=7            
    
    I_Kind=8            !IJKBC(:,I_Kind),�������2�е�8�ж���Ϊ�߽���������
    I_Kindsub=9         !IJKBC(:,I_Kindsub),�������2�е�9�ж��嶨��߽���������
    
    I_Blk_Ex=10         !IJKBC(:,I_Blk_Ex),�������2�е�10�ж��嶨���Ե��ӶԷ����      
    I_Blkbd_Ex=11       !IJKBC(:,I_Blkbd_Ex),�������2�е�11�ж��嶨���Ե��ӶԷ���߽��
    I_Lcross=12         !IJKBC(:,I_Lcross),�������2�е�12�ж��嶨�屾�������Ƿ����Ӧ�Է����Ƿ񽻲棨����������߽�����I�����Ӧ�Է������߽����ϵ�I������J����
    I_Lrev_A=13         !IJKBC(:,I_Lrev_A),�������2�е�13�ж��嶨�屾������Ƿ����Ӧ�Է�����I���������棨���������߽�����I��Сֵ��Ӧ�Է������I���ֵ������Сֵ��
    I_Lrev_B=14         !IJKBC(:,I_Lrev_B),�������2�е�14�ж��嶨�屾������Ƿ����Ӧ�Է�����J���������棨���������߽�����I��Сֵ��Ӧ�Է������J���ֵ������Сֵ��
   I_Period=15

END SUBROUTINE COMMON_Define_BCinput_Parameter