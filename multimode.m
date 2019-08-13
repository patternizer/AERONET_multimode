function [GMM]=multimode(Y)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SYNTAX:                                                       %
    % function [GMM]=multimode(Y)                                   %
    %                                                               %
    % REFERENCE:                                                    %
    % Taylor, M., Kazadzis, S., Gerosopoulos, E. (2014)             % 
    % "Multi-modal analysis of aerosol robotic network              % 
    % size distributions for remote sensing applications"           % 
    % Atmospheric Measurement Techniques 7, 839-858.                % 
    % doi:10.5194/amt-7-839-2014                                    %
    %                                                               %
    % INPUTS:                                                       %
    % Y  --> 22-bin AERONET aerosol size distribution function      %
    %                                                               %
    % OUTPUT:                                                       %
    % GMM --> structure array containing:                           %
    % GMM.stats --> mode,b,s,R2,V                                   %
    % GMM.modes --> modal V(1-8), r(1-8), sigma(1-8)                %
    % GMM.coef  --> lognormal coefficients: a(1-8), b(1-8), c(1-8)  %
    % GMM.se    --> standard errors on V,r,sigma (1-8)              %
    % GMM.test  --> Fisher test statistics                          %
    %                                                               %
    % DEPENDENCIES:                                                 %
    % none                                                          %
    %                                                               %
    % Version 1.0 (06/05/2016)                                      %
    % Dr Michael Taylor                                             %
    % LAP/AUTh: http://lap.physics.auth.gr/faculty.asp?lang=en      %
    % URL: http://users.auth.gr/mtaylor/                            %
    % Email (1): mtaylor@auth.gr                                    %
    % Email (2): patternizer@gmail.com                              %
    % ORCID: http://orcid.org/0000-0002-3473-3478                   %
    % ResearchGate: http://bit.ly/1WWrLoJ                           %
    % GoogleScholar: http://bit.ly/1LAK8io                          %
    %                                                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 	% clear all; close all; clc;
   
    %% I/O PARAMETERS
    path_save           = ''; 
    
    %% PLOT PARAMETERS
    line_width          = 2;
    marker_size         = 10;
    font_size           = 18;
    font_weight         = 'bold';           
    grey                = [204 204 204]/255;    % Light Graphite
    bright              = [255 0 0]/255;        % Red
    mid                 = [237 138 237]/255;    % Violet
    dark                = [138 0 138]/255;      % Dark Magenta
    opposite1           = [0 0 255]/255;        % Blue
    opposite2           = [0 0 138]/255;        % Navy 
    opposite3           = [0 204 204]/255;      % Green
    figure_type         = 2;                    % 1=FIG, 2=PNG     
    figure_count        = 0;    
    FLAG_EPS            = 0;                    % 0=none, 1=generate EPS
    FLAG_plot_on        = 1;    
    FLAG_plot_visible   = 1;    
    FLAG_plot_save      = 1;    
    FLAG_plot_GMM       = 1;

    %% MODEL PARAMETERS
    FLAG_autodetect     = 1; % 0=specify uppper limit to n-modes, 1=auto (1:8 modes)
    FLAG_method         = 2; % 0=Multiple Regression, 1=Fisher, 2=Subjective
    if isempty(Y)    
        % AERONET Aerosol Inversions L1.5 (Version 2): EXETER_MO, 28/01/2016, 10:29:08 (GMT)    
        Y = [0.000262,0.001774,0.005377,0.007652,0.005960,0.003439,0.002175,0.002049,0.003043,0.005746,0.009857,0.013233,0.016114,0.019276,0.022546,0.024126,0.022326,0.016942,0.009724,0.003899,0.001043,0.000184];
    else
    end    
    X                   = [ 0.050000 , 0.065604 , 0.086077 , 0.112939 , 0.148184 , 0.194429 , 0.255105 , 0.334716 , 0.439173 , 0.576227 , 0.756052 , 0.991996 , 1.301571 , 1.707757 , 2.240702 , 2.939966 , 3.857452 , 5.061260 , 6.640745 , 8.713145 , 11.432287 , 15.00000 ];
    str_x               = 'r [\mum]';
    str_y               = 'dV/dlnr [\mum^3/\mum^2]';    
    str_Y               = 'V';
    str_precision       = '%5.3f'; 
    N_interp            = int32(22*(2^3)); % interpolation doubling factor
    n_interp            = double(N_interp);
    
    %% SMOOTH PIECEWISE INTERPOLATION [n_interp]
    ft = fittype('pchipinterp'); 
    opts = fitoptions(ft);
    [fitresult,gof] = fit(X',Y',ft,opts);
    xmin = log10(min(X)); 
    xmax = log10(max(X)); 
    x_vector = logspace(xmin,xmax,N_interp);
    y_vector = fitresult(x_vector)';
    ymax = max(y_vector);
    V_trapz = trapz(log(X),Y);
    V_interp = trapz(log(x_vector),y_vector);
    X_interp = x_vector'; % r
    Y_interp = y_vector'; % dV/dln(r)
    logX = log(x_vector)'; % ln(r)
    clear x_vector y_vector
        
    %% MULTI-MODAL GUASSIAN FITS 1-8

    MATRIX_GMM          = ones(length(Y),31)*NaN;
    GMM_coef    = zeros(8,25);
    GMM_se       = zeros(8,25);

        a_min = 1E-3;
        b_min = log(0.05);
        c_min = sqrt(2)*0.1;
        
        % GAUSS-1 FIT
        ft = fittype('gauss1');
        opts = fitoptions(ft);
        a1_min=a_min; b1_min=b_min; c1_min=c_min;        
        a1_max=inf; b1_max=inf; c1_max=inf;
        opts.Lower = [a1_min b1_min c1_min];
        opts.Upper = [a1_max b1_max c1_max];
        opts.Robust = 'Bisquare';
        [fitresult1,gof1] = fit(logX,Y_interp,ft,opts);
        y_G1 = fitresult1(logX);            
        A1_opt_G1=fitresult1.a1;
        B1_opt_G1=fitresult1.b1;
        C1_opt_G1=fitresult1.c1;        
        df_a1=exp(-((logX-B1_opt_G1).^2./C1_opt_G1^2));
        df_b1=exp(-((logX-B1_opt_G1).^2./C1_opt_G1^2)) .* 2.*A1_opt_G1.*(logX-B1_opt_G1)./C1_opt_G1^2;
        df_c1=exp(-((logX-B1_opt_G1).^2./C1_opt_G1^2)) .* 2.*A1_opt_G1.*(logX-B1_opt_G1).^2./C1_opt_G1^3;                      
        ci = confint(fitresult1,0.95); ci(isnan(ci))=0;
        SE_a1=(A1_opt_G1-ci(1,1))/1.96;
        SE_b1=(B1_opt_G1-ci(1,2))/1.96;
        SE_c1=(C1_opt_G1-ci(1,3))/1.96;
        SE_G1=abs(sqrt( (df_a1.*SE_a1).^2+(df_b1.*SE_b1).^2+(df_c1.*SE_c1).^2 ));
        y_G1_U=y_G1+SE_G1;
        y_G1_L=y_G1-SE_G1;
        V_G1 = trapz(logX,y_G1); str_G1={['1-modal fit: ',str_Y,'=',num2str(V_G1,str_precision)]};
        y_G1_U=y_G1+SE_G1; y_G1_L=y_G1-SE_G1;        
        GMM_coef(1,:)=[1,A1_opt_G1,zeros(1,7),B1_opt_G1,zeros(1,7),C1_opt_G1,zeros(1,7)];
        GMM_se(1,:)=[1,SE_a1,zeros(1,7),SE_b1,zeros(1,7),SE_c1,zeros(1,7)];

        % GAUSS-2 FIT
        ft = fittype('gauss2');
        opts = fitoptions(ft);  
        a1_min=a_min; b1_min=b_min; c1_min=c_min;                
        a2_min=a_min; b2_min=b_min; c2_min=c_min;                
        a1_max=inf; b1_max=inf; c1_max=inf;                
        a2_max=inf; b2_max=inf; c2_max=inf;
        opts.Lower = [a1_min b1_min c1_min a2_min b2_min c2_min];        
        opts.Upper = [a1_max b1_max c1_max a2_max b2_max c2_max];        
        opts.Robust = 'Bisquare';
        [fitresult2,gof2] = fit(logX,Y_interp,ft,opts);
        y_G2 = fitresult2(logX);
        A1_opt_G2=fitresult2.a1;
        A2_opt_G2=fitresult2.a2;
        B1_opt_G2=fitresult2.b1;
        B2_opt_G2=fitresult2.b2;
        C1_opt_G2=fitresult2.c1;
        C2_opt_G2=fitresult2.c2;        
        df_a1=exp(-((logX-B1_opt_G2).^2./C1_opt_G2^2));
        df_a2=exp(-((logX-B2_opt_G2).^2./C2_opt_G2^2));
        df_b1=exp(-((logX-B1_opt_G2).^2./C1_opt_G2^2)) .* 2.*A1_opt_G2.*(logX-B1_opt_G2)./C1_opt_G2^2;
        df_b2=exp(-((logX-B2_opt_G2).^2./C2_opt_G2^2)) .* 2.*A2_opt_G2.*(logX-B2_opt_G2)./C2_opt_G2^2;
        df_c1=exp(-((logX-B1_opt_G2).^2./C1_opt_G2^2)) .* 2.*A1_opt_G2.*(logX-B1_opt_G2).^2./C1_opt_G2^3;                      
        df_c2=exp(-((logX-B2_opt_G2).^2./C2_opt_G2^2)) .* 2.*A2_opt_G2.*(logX-B2_opt_G2).^2./C2_opt_G2^3;                      
        ci = confint(fitresult2,0.95); ci(isnan(ci))=0;
        SE_a1=(A1_opt_G2-ci(1,1))/1.96;
        SE_b1=(B1_opt_G2-ci(1,2))/1.96;
        SE_c1=(C1_opt_G2-ci(1,3))/1.96;
        SE_a2=(A2_opt_G2-ci(1,4))/1.96;
        SE_b2=(B2_opt_G2-ci(1,5))/1.96;
        SE_c2=(C2_opt_G2-ci(1,6))/1.96;
        SE_G2=abs(sqrt( (df_a1.*SE_a1).^2+(df_b1.*SE_b1).^2+(df_c1.*SE_c1).^2 + (df_a2.*SE_a2).^2+(df_b2.*SE_b2).^2+(df_c2.*SE_c2).^2 ));
        y_G2_U=y_G2+SE_G2;
        y_G2_L=y_G2-SE_G2;        
        V_G2 = trapz(logX,y_G2); str_G2={['2-modal fit: ',str_Y,'=',num2str(V_G2,str_precision)]}; 
        y_G2_U=y_G2+SE_G2; y_G2_L=y_G2-SE_G2;
        GMM_coef(2,:)=[2,A1_opt_G2,A2_opt_G2,zeros(1,6),B1_opt_G2,B2_opt_G2,zeros(1,6),C1_opt_G2,C2_opt_G2,zeros(1,6)];
        GMM_se(2,:)=[2,SE_a1,SE_a2,zeros(1,6),SE_b1,SE_b2,zeros(1,6),SE_c1,SE_c2,zeros(1,6)];

        % GAUSS-3 FIT
        ft = fittype( 'gauss3' );
        opts = fitoptions( ft );             
        a1_min=a_min; b1_min=b_min; c1_min=c_min;                
        a2_min=a_min; b2_min=b_min; c2_min=c_min;               
        a3_min=a_min; b3_min=b_min; c3_min=c_min;                
        a1_max=inf; b1_max=inf; c1_max=inf;                
        a2_max=inf; b2_max=inf; c2_max=inf;                
        a3_max=inf; b3_max=inf; c3_max=inf;  
        opts.Lower = [a1_min b1_min c1_min a2_min b2_min c2_min a3_min b3_min c3_min];
        opts.Upper = [a1_max b1_max c1_max a2_max b2_max c2_max a3_max b3_max c3_max];
        opts.Robust = 'Bisquare';
        [fitresult3,gof3] = fit(logX,Y_interp,ft,opts);
        y_G3 = fitresult3(logX);
        A1_opt_G3=fitresult3.a1;
        A2_opt_G3=fitresult3.a2;
        A3_opt_G3=fitresult3.a3;
        B1_opt_G3=fitresult3.b1;
        B2_opt_G3=fitresult3.b2;
        B3_opt_G3=fitresult3.b3;
        C1_opt_G3=fitresult3.c1;
        C2_opt_G3=fitresult3.c2;
        C3_opt_G3=fitresult3.c3;
        df_a1=exp(-((logX-B1_opt_G3).^2./C1_opt_G3^2));
        df_a2=exp(-((logX-B2_opt_G3).^2./C2_opt_G3^2));
        df_a3=exp(-((logX-B3_opt_G3).^2./C3_opt_G3^2));
        df_b1=exp(-((logX-B1_opt_G3).^2./C1_opt_G3^2)) .* 2.*A1_opt_G3.*(logX-B1_opt_G3)./C1_opt_G3^2;
        df_b2=exp(-((logX-B2_opt_G3).^2./C2_opt_G3^2)) .* 2.*A2_opt_G3.*(logX-B2_opt_G3)./C2_opt_G3^2;
        df_b3=exp(-((logX-B3_opt_G3).^2./C3_opt_G3^2)) .* 2.*A3_opt_G3.*(logX-B3_opt_G3)./C3_opt_G3^2;
        df_c1=exp(-((logX-B1_opt_G3).^2./C1_opt_G3^2)) .* 2.*A1_opt_G3.*(logX-B1_opt_G3).^2./C1_opt_G3^3;                      
        df_c2=exp(-((logX-B2_opt_G3).^2./C2_opt_G3^2)) .* 2.*A2_opt_G3.*(logX-B2_opt_G3).^2./C2_opt_G3^3;                      
        df_c3=exp(-((logX-B3_opt_G3).^2./C3_opt_G3^2)) .* 2.*A3_opt_G3.*(logX-B3_opt_G3).^2./C3_opt_G3^3;                      
        ci = confint(fitresult3,0.95); ci(isnan(ci))=0;
        SE_a1=(A1_opt_G3-ci(1,1))/1.96;
        SE_b1=(B1_opt_G3-ci(1,2))/1.96;
        SE_c1=(C1_opt_G3-ci(1,3))/1.96;
        SE_a2=(A2_opt_G3-ci(1,4))/1.96;
        SE_b2=(B2_opt_G3-ci(1,5))/1.96;
        SE_c2=(C2_opt_G3-ci(1,6))/1.96;
        SE_a3=(A3_opt_G3-ci(1,7))/1.96;
        SE_b3=(B3_opt_G3-ci(1,8))/1.96;
        SE_c3=(C3_opt_G3-ci(1,9))/1.96;
        SE_G3=abs(sqrt( (df_a1.*SE_a1).^2+(df_b1.*SE_b1).^2+(df_c1.*SE_c1).^2 + (df_a2.*SE_a2).^2+(df_b2.*SE_b2).^2+(df_c2.*SE_c2).^2 + (df_a3.*SE_a3).^2+(df_b3.*SE_b3).^2+(df_c3.*SE_c3).^2    ));
        y_G3_U=y_G3+SE_G3;
        y_G3_L=y_G3-SE_G3;  
        V_G3 = trapz(logX,y_G3); str_G3={['3-modal fit: ',str_Y,'=',num2str(V_G3,str_precision)]}; 
        y_G3_U=y_G3+SE_G3; y_G3_L=y_G3-SE_G3;
        GMM_coef(3,:)=[3,A1_opt_G3,A2_opt_G3,A3_opt_G3,zeros(1,5),B1_opt_G3,B2_opt_G3,B3_opt_G3,zeros(1,5),C1_opt_G3,C2_opt_G3,C3_opt_G3,zeros(1,5)];
        GMM_se(3,:)=[3,SE_a1,SE_a2,SE_a3,zeros(1,5),SE_b1,SE_b2,SE_b3,zeros(1,5),SE_c1,SE_c2,SE_c3,zeros(1,5)];
      
        % GAUSS-4 FIT
        ft = fittype( 'gauss4' );
        opts = fitoptions( ft );             
        a1_min=a_min; b1_min=b_min; c1_min=c_min;                
        a2_min=a_min; b2_min=b_min; c2_min=c_min;                
        a3_min=a_min; b3_min=b_min; c3_min=c_min;               
        a4_min=a_min; b4_min=b_min; c4_min=c_min;                
        a1_max=inf; b1_max=inf; c1_max=inf;                
        a2_max=inf; b2_max=inf; c2_max=inf;                
        a3_max=inf; b3_max=inf; c3_max=inf;          
        a4_max=inf; b4_max=inf; c4_max=inf;  
        opts.Lower = [a1_min b1_min c1_min a2_min b2_min c2_min a3_min b3_min c3_min a4_min b4_min c4_min];
        opts.Upper = [a1_max b1_max c1_max a2_max b2_max c2_max a3_max b3_max c3_max a4_max b4_max c4_max];
        opts.Robust = 'Bisquare';
        [fitresult4,gof4] = fit(logX,Y_interp,ft,opts);
        y_G4 = fitresult4(logX);
        A1_opt_G4=fitresult4.a1;
        A2_opt_G4=fitresult4.a2;
        A3_opt_G4=fitresult4.a3;
        A4_opt_G4=fitresult4.a4;
        B1_opt_G4=fitresult4.b1;
        B2_opt_G4=fitresult4.b2;
        B3_opt_G4=fitresult4.b3;
        B4_opt_G4=fitresult4.b4;
        C1_opt_G4=fitresult4.c1;
        C2_opt_G4=fitresult4.c2;
        C3_opt_G4=fitresult4.c3;    
        C4_opt_G4=fitresult4.c4;    
        df_a1=exp(-((logX-B1_opt_G4).^2./C1_opt_G4^2));
        df_a2=exp(-((logX-B2_opt_G4).^2./C2_opt_G4^2));
        df_a3=exp(-((logX-B3_opt_G4).^2./C3_opt_G4^2));
        df_a4=exp(-((logX-B4_opt_G4).^2./C4_opt_G4^2));
        df_b1=exp(-((logX-B1_opt_G4).^2./C1_opt_G4^2)) .* 2.*A1_opt_G4.*(logX-B1_opt_G4)./C1_opt_G4^2;
        df_b2=exp(-((logX-B2_opt_G4).^2./C2_opt_G4^2)) .* 2.*A2_opt_G4.*(logX-B2_opt_G4)./C2_opt_G4^2;
        df_b3=exp(-((logX-B3_opt_G4).^2./C3_opt_G4^2)) .* 2.*A3_opt_G4.*(logX-B3_opt_G4)./C3_opt_G4^2;
        df_b4=exp(-((logX-B4_opt_G4).^2./C4_opt_G4^2)) .* 2.*A4_opt_G4.*(logX-B4_opt_G4)./C4_opt_G4^2;
        df_c1=exp(-((logX-B1_opt_G4).^2./C1_opt_G4^2)) .* 2.*A1_opt_G4.*(logX-B1_opt_G4).^2./C1_opt_G4^3;                      
        df_c2=exp(-((logX-B2_opt_G4).^2./C2_opt_G4^2)) .* 2.*A2_opt_G4.*(logX-B2_opt_G4).^2./C2_opt_G4^3;                      
        df_c3=exp(-((logX-B3_opt_G4).^2./C3_opt_G4^2)) .* 2.*A3_opt_G4.*(logX-B3_opt_G4).^2./C3_opt_G4^3;                      
        df_c4=exp(-((logX-B4_opt_G4).^2./C4_opt_G4^2)) .* 2.*A4_opt_G4.*(logX-B4_opt_G4).^2./C4_opt_G4^3;                      
        ci = confint(fitresult4,0.95); ci(isnan(ci))=0;
        SE_a1=(A1_opt_G4-ci(1,1))/1.96;
        SE_b1=(B1_opt_G4-ci(1,2))/1.96;
        SE_c1=(C1_opt_G4-ci(1,3))/1.96;
        SE_a2=(A2_opt_G4-ci(1,4))/1.96;
        SE_b2=(B2_opt_G4-ci(1,5))/1.96;
        SE_c2=(C2_opt_G4-ci(1,6))/1.96;
        SE_a3=(A3_opt_G4-ci(1,7))/1.96;
        SE_b3=(B3_opt_G4-ci(1,8))/1.96;
        SE_c3=(C3_opt_G4-ci(1,9))/1.96;
        SE_a4=(A4_opt_G4-ci(1,10))/1.96;
        SE_b4=(B4_opt_G4-ci(1,11))/1.96;
        SE_c4=(C4_opt_G4-ci(1,12))/1.96;
        SE_G4=abs(sqrt( (df_a1.*SE_a1).^2+(df_b1.*SE_b1).^2+(df_c1.*SE_c1).^2 + (df_a2.*SE_a2).^2+(df_b2.*SE_b2).^2+(df_c2.*SE_c2).^2 + (df_a3.*SE_a3).^2+(df_b3.*SE_b3).^2+(df_c3.*SE_c3).^2  + (df_a4.*SE_a4).^2+(df_b4.*SE_b4).^2+(df_c4.*SE_c4).^2    ));
        y_G4_U=y_G4+SE_G4;
        y_G4_L=y_G4-SE_G4; 
        V_G4 = trapz(logX,y_G4); str_G4={['4-modal fit: ',str_Y,'=',num2str(V_G4,str_precision)]};
        y_G4_U=y_G4+SE_G4; y_G4_L=y_G4-SE_G4;
        GMM_coef(4,:)=[4,A1_opt_G4,A2_opt_G4,A3_opt_G4,A4_opt_G4,zeros(1,4),B1_opt_G4,B2_opt_G4,B3_opt_G4,B4_opt_G4,zeros(1,4),C1_opt_G4,C2_opt_G4,C3_opt_G4,C4_opt_G4,zeros(1,4)];
        GMM_se(4,:)=[4,SE_a1,SE_a2,SE_a3,SE_a4,zeros(1,4),SE_b1,SE_b2,SE_b3,SE_b4,zeros(1,4),SE_c1,SE_c2,SE_c3,SE_c4,zeros(1,4)];
        
        % GAUSS-5 FIT
        ft = fittype( 'gauss5' );
        opts = fitoptions( ft );            
        a1_min=a_min; b1_min=b_min; c1_min=c_min;        
        a2_min=a_min; b2_min=b_min; c2_min=c_min;                
        a3_min=a_min; b3_min=b_min; c3_min=c_min;                
        a4_min=a_min; b4_min=b_min; c4_min=c_min;                
        a5_min=a_min; b5_min=b_min; c5_min=c_min;                
        a1_max=inf; b1_max=inf; c1_max=inf;                
        a2_max=inf; b2_max=inf; c2_max=inf;                
        a3_max=inf; b3_max=inf; c3_max=inf;          
        a4_max=inf; b4_max=inf; c4_max=inf;         
        a5_max=inf; b5_max=inf; c5_max=inf; 
        opts.Lower = [a1_min b1_min c1_min a2_min b2_min c2_min a3_min b3_min c3_min a4_min b4_min c4_min a5_min b5_min c5_min];
        opts.Upper = [a1_max b1_max c1_max a2_max b2_max c2_max a3_max b3_max c3_max a4_max b4_max c4_max a5_max b5_max c5_max];
        opts.Robust = 'Bisquare';
        [fitresult5,gof5] = fit(logX,Y_interp,ft,opts);
        y_G5 = fitresult5(logX);
        A1_opt_G5=fitresult5.a1;
        A2_opt_G5=fitresult5.a2;
        A3_opt_G5=fitresult5.a3;
        A4_opt_G5=fitresult5.a4;
        A5_opt_G5=fitresult5.a5;
        B1_opt_G5=fitresult5.b1;
        B2_opt_G5=fitresult5.b2;
        B3_opt_G5=fitresult5.b3;
        B4_opt_G5=fitresult5.b4;
        B5_opt_G5=fitresult5.b5;
        C1_opt_G5=fitresult5.c1;
        C2_opt_G5=fitresult5.c2;
        C3_opt_G5=fitresult5.c3;    
        C4_opt_G5=fitresult5.c4;    
        C5_opt_G5=fitresult5.c5;         
        df_a1=exp(-((logX-B1_opt_G5).^2./C1_opt_G5^2));
        df_a2=exp(-((logX-B2_opt_G5).^2./C2_opt_G5^2));
        df_a3=exp(-((logX-B3_opt_G5).^2./C3_opt_G5^2));
        df_a4=exp(-((logX-B4_opt_G5).^2./C4_opt_G5^2));
        df_a5=exp(-((logX-B5_opt_G5).^2./C5_opt_G5^2));
        df_b1=exp(-((logX-B1_opt_G5).^2./C1_opt_G5^2)) .* 2.*A1_opt_G5.*(logX-B1_opt_G5)./C1_opt_G5^2;
        df_b2=exp(-((logX-B2_opt_G5).^2./C2_opt_G5^2)) .* 2.*A2_opt_G5.*(logX-B2_opt_G5)./C2_opt_G5^2;
        df_b3=exp(-((logX-B3_opt_G5).^2./C3_opt_G5^2)) .* 2.*A3_opt_G5.*(logX-B3_opt_G5)./C3_opt_G5^2;
        df_b4=exp(-((logX-B4_opt_G5).^2./C4_opt_G5^2)) .* 2.*A4_opt_G5.*(logX-B4_opt_G5)./C4_opt_G5^2;
        df_b5=exp(-((logX-B5_opt_G5).^2./C5_opt_G5^2)) .* 2.*A5_opt_G5.*(logX-B5_opt_G5)./C5_opt_G5^2;
        df_c1=exp(-((logX-B1_opt_G5).^2./C1_opt_G5^2)) .* 2.*A1_opt_G5.*(logX-B1_opt_G5).^2./C1_opt_G5^3;                      
        df_c2=exp(-((logX-B2_opt_G5).^2./C2_opt_G5^2)) .* 2.*A2_opt_G5.*(logX-B2_opt_G5).^2./C2_opt_G5^3;                      
        df_c3=exp(-((logX-B3_opt_G5).^2./C3_opt_G5^2)) .* 2.*A3_opt_G5.*(logX-B3_opt_G5).^2./C3_opt_G5^3;                      
        df_c4=exp(-((logX-B4_opt_G5).^2./C4_opt_G5^2)) .* 2.*A4_opt_G5.*(logX-B4_opt_G5).^2./C4_opt_G5^3;                      
        df_c5=exp(-((logX-B5_opt_G5).^2./C5_opt_G5^2)) .* 2.*A5_opt_G5.*(logX-B5_opt_G5).^2./C5_opt_G5^3;                      
        ci = confint(fitresult5,0.95); ci(isnan(ci))=0;
        SE_a1=(A1_opt_G5-ci(1,1))/1.96;
        SE_b1=(B1_opt_G5-ci(1,2))/1.96;
        SE_c1=(C1_opt_G5-ci(1,3))/1.96;
        SE_a2=(A2_opt_G5-ci(1,4))/1.96;
        SE_b2=(B2_opt_G5-ci(1,5))/1.96;
        SE_c2=(C2_opt_G5-ci(1,6))/1.96;
        SE_a3=(A3_opt_G5-ci(1,7))/1.96;
        SE_b3=(B3_opt_G5-ci(1,8))/1.96;
        SE_c3=(C3_opt_G5-ci(1,9))/1.96;
        SE_a4=(A4_opt_G5-ci(1,10))/1.96;
        SE_b4=(B4_opt_G5-ci(1,11))/1.96;
        SE_c4=(C4_opt_G5-ci(1,12))/1.96;
        SE_a5=(A5_opt_G5-ci(1,13))/1.96;
        SE_b5=(B5_opt_G5-ci(1,14))/1.96;
        SE_c5=(C5_opt_G5-ci(1,15))/1.96;
        SE_G5=abs(sqrt( (df_a1.*SE_a1).^2+(df_b1.*SE_b1).^2+(df_c1.*SE_c1).^2 + (df_a2.*SE_a2).^2+(df_b2.*SE_b2).^2+(df_c2.*SE_c2).^2 + (df_a3.*SE_a3).^2+(df_b3.*SE_b3).^2+(df_c3.*SE_c3).^2  + (df_a4.*SE_a4).^2+(df_b4.*SE_b4).^2+(df_c4.*SE_c4).^2  + (df_a5.*SE_a5).^2+(df_b5.*SE_b5).^2+(df_c5.*SE_c5).^2    ));
        y_G5_U=y_G5+SE_G5;
        y_G5_L=y_G5-SE_G5;         
        V_G5 = trapz(logX,y_G5); str_G5={['5-modal fit: ',str_Y,'=',num2str(V_G5,str_precision)]}; 
        y_G5_U=y_G5+SE_G5; y_G5_L=y_G5-SE_G5;
        GMM_coef(5,:)=[5,A1_opt_G5,A2_opt_G5,A3_opt_G5,A4_opt_G5,A5_opt_G5,zeros(1,3),B1_opt_G5,B2_opt_G5,B3_opt_G5,B4_opt_G5,B5_opt_G5,zeros(1,3),C1_opt_G5,C2_opt_G5,C3_opt_G5,C4_opt_G5,C5_opt_G5,zeros(1,3)];
        GMM_se(5,:)=[5,SE_a1,SE_a2,SE_a3,SE_a4,SE_a5,zeros(1,3),SE_b1,SE_b2,SE_b3,SE_b4,SE_b5,zeros(1,3),SE_c1,SE_c2,SE_c3,SE_c4,SE_c5,zeros(1,3)];
 
        % GAUSS-6 FIT
        ft = fittype( 'gauss6' );
        opts = fitoptions( ft );           
        a1_min=a_min; b1_min=b_min; c1_min=c_min;                
        a2_min=a_min; b2_min=b_min; c2_min=c_min;                
        a3_min=a_min; b3_min=b_min; c3_min=c_min;                
        a4_min=a_min; b4_min=b_min; c4_min=c_min;                
        a5_min=a_min; b5_min=b_min; c5_min=c_min;                
        a6_min=a_min; b6_min=b_min; c6_min=c_min;                
        a1_max=inf; b1_max=inf; c1_max=inf;                
        a2_max=inf; b2_max=inf; c2_max=inf;                
        a3_max=inf; b3_max=inf; c3_max=inf;          
        a4_max=inf; b4_max=inf; c4_max=inf;         
        a5_max=inf; b5_max=inf; c5_max=inf;        
        a6_max=inf; b6_max=inf; c6_max=inf; 
        opts.Lower = [a1_min b1_min c1_min a2_min b2_min c2_min a3_min b3_min c3_min a4_min b4_min c4_min a5_min b5_min c5_min a6_min b6_min c6_min];
        opts.Upper = [a1_max b1_max c1_max a2_max b2_max c2_max a3_max b3_max c3_max a4_max b4_max c4_max a5_max b5_max c5_max a6_max b6_max c6_max];
        opts.Robust = 'Bisquare';
        [fitresult6,gof6] = fit(logX,Y_interp,ft,opts);
        y_G6 = fitresult6(logX);
        A1_opt_G6=fitresult6.a1;
        A2_opt_G6=fitresult6.a2;
        A3_opt_G6=fitresult6.a3;
        A4_opt_G6=fitresult6.a4;
        A5_opt_G6=fitresult6.a5;
        A6_opt_G6=fitresult6.a6;
        B1_opt_G6=fitresult6.b1;
        B2_opt_G6=fitresult6.b2;
        B3_opt_G6=fitresult6.b3;
        B4_opt_G6=fitresult6.b4;
        B5_opt_G6=fitresult6.b5;
        B6_opt_G6=fitresult6.b6;
        C1_opt_G6=fitresult6.c1;
        C2_opt_G6=fitresult6.c2;
        C3_opt_G6=fitresult6.c3;    
        C4_opt_G6=fitresult6.c4;    
        C5_opt_G6=fitresult6.c5;    
        C6_opt_G6=fitresult6.c6;                            
        df_a1=exp(-((logX-B1_opt_G6).^2./C1_opt_G6^2));
        df_a2=exp(-((logX-B2_opt_G6).^2./C2_opt_G6^2));
        df_a3=exp(-((logX-B3_opt_G6).^2./C3_opt_G6^2));
        df_a4=exp(-((logX-B4_opt_G6).^2./C4_opt_G6^2));
        df_a5=exp(-((logX-B5_opt_G6).^2./C5_opt_G6^2));
        df_a6=exp(-((logX-B6_opt_G6).^2./C6_opt_G6^2));
        df_b1=exp(-((logX-B1_opt_G6).^2./C1_opt_G6^2)) .* 2.*A1_opt_G6.*(logX-B1_opt_G6)./C1_opt_G6^2;
        df_b2=exp(-((logX-B2_opt_G6).^2./C2_opt_G6^2)) .* 2.*A2_opt_G6.*(logX-B2_opt_G6)./C2_opt_G6^2;
        df_b3=exp(-((logX-B3_opt_G6).^2./C3_opt_G6^2)) .* 2.*A3_opt_G6.*(logX-B3_opt_G6)./C3_opt_G6^2;
        df_b4=exp(-((logX-B4_opt_G6).^2./C4_opt_G6^2)) .* 2.*A4_opt_G6.*(logX-B4_opt_G6)./C4_opt_G6^2;
        df_b5=exp(-((logX-B5_opt_G6).^2./C5_opt_G6^2)) .* 2.*A5_opt_G6.*(logX-B5_opt_G6)./C5_opt_G6^2;
        df_b6=exp(-((logX-B6_opt_G6).^2./C6_opt_G6^2)) .* 2.*A6_opt_G6.*(logX-B6_opt_G6)./C6_opt_G6^2;
        df_c1=exp(-((logX-B1_opt_G6).^2./C1_opt_G6^2)) .* 2.*A1_opt_G6.*(logX-B1_opt_G6).^2./C1_opt_G6^3;                      
        df_c2=exp(-((logX-B2_opt_G6).^2./C2_opt_G6^2)) .* 2.*A2_opt_G6.*(logX-B2_opt_G6).^2./C2_opt_G6^3;                      
        df_c3=exp(-((logX-B3_opt_G6).^2./C3_opt_G6^2)) .* 2.*A3_opt_G6.*(logX-B3_opt_G6).^2./C3_opt_G6^3;                      
        df_c4=exp(-((logX-B4_opt_G6).^2./C4_opt_G6^2)) .* 2.*A4_opt_G6.*(logX-B4_opt_G6).^2./C4_opt_G6^3;                      
        df_c5=exp(-((logX-B5_opt_G6).^2./C5_opt_G6^2)) .* 2.*A5_opt_G6.*(logX-B5_opt_G6).^2./C5_opt_G6^3;                      
        df_c6=exp(-((logX-B6_opt_G6).^2./C6_opt_G6^2)) .* 2.*A6_opt_G6.*(logX-B6_opt_G6).^2./C6_opt_G6^3;                      
        ci = confint(fitresult6,0.95); ci(isnan(ci))=0;
        SE_a1=(A1_opt_G6-ci(1,1))/1.96;
        SE_b1=(B1_opt_G6-ci(1,2))/1.96;
        SE_c1=(C1_opt_G6-ci(1,3))/1.96;
        SE_a2=(A2_opt_G6-ci(1,4))/1.96;
        SE_b2=(B2_opt_G6-ci(1,5))/1.96;
        SE_c2=(C2_opt_G6-ci(1,6))/1.96;
        SE_a3=(A3_opt_G6-ci(1,7))/1.96;
        SE_b3=(B3_opt_G6-ci(1,8))/1.96;
        SE_c3=(C3_opt_G6-ci(1,9))/1.96;
        SE_a4=(A4_opt_G6-ci(1,10))/1.96;
        SE_b4=(B4_opt_G6-ci(1,11))/1.96;
        SE_c4=(C4_opt_G6-ci(1,12))/1.96;
        SE_a5=(A5_opt_G6-ci(1,13))/1.96;
        SE_b5=(B5_opt_G6-ci(1,14))/1.96;
        SE_c5=(C5_opt_G6-ci(1,15))/1.96;
        SE_a6=(A6_opt_G6-ci(1,16))/1.96;
        SE_b6=(B6_opt_G6-ci(1,17))/1.96;
        SE_c6=(C6_opt_G6-ci(1,18))/1.96;        
        SE_G6=abs(sqrt( (df_a1.*SE_a1).^2+(df_b1.*SE_b1).^2+(df_c1.*SE_c1).^2 + (df_a2.*SE_a2).^2+(df_b2.*SE_b2).^2+(df_c2.*SE_c2).^2 + (df_a3.*SE_a3).^2+(df_b3.*SE_b3).^2+(df_c3.*SE_c3).^2  + (df_a4.*SE_a4).^2+(df_b4.*SE_b4).^2+(df_c4.*SE_c4).^2  + (df_a5.*SE_a5).^2+(df_b5.*SE_b5).^2+(df_c5.*SE_c5).^2  + (df_a6.*SE_a6).^2+(df_b6.*SE_b6).^2+(df_c6.*SE_c6).^2    ));
        y_G6_U=y_G6+SE_G6;
        y_G6_L=y_G6-SE_G6;   
        V_G6 = trapz(logX,y_G6); str_G6={['6-modal fit: ',str_Y,'=',num2str(V_G6,str_precision)]};
        y_G6_U=y_G6+SE_G6; y_G6_L=y_G6-SE_G6;
        GMM_coef(6,:)=[6,A1_opt_G6,A2_opt_G6,A3_opt_G6,A4_opt_G6,A5_opt_G6,A6_opt_G6,zeros(1,2),B1_opt_G6,B2_opt_G6,B3_opt_G6,B4_opt_G6,B5_opt_G6,B6_opt_G6,zeros(1,2),C1_opt_G6,C2_opt_G6,C3_opt_G6,C4_opt_G6,C5_opt_G6,C6_opt_G6,zeros(1,2)];
        GMM_se(6,:)=[6,SE_a1,SE_a2,SE_a3,SE_a4,SE_a5,SE_a6,zeros(1,2),SE_b1,SE_b2,SE_b3,SE_b4,SE_b5,SE_b6,zeros(1,2),SE_c1,SE_c2,SE_c3,SE_c4,SE_c5,SE_c6,zeros(1,2)];
       
        % GAUSS-7 FIT
        ft = fittype( 'gauss7' );
        opts = fitoptions( ft );             
        a1_min=a_min; b1_min=b_min; c1_min=c_min;                
        a2_min=a_min; b2_min=b_min; c2_min=c_min;                
        a3_min=a_min; b3_min=b_min; c3_min=c_min;                
        a4_min=a_min; b4_min=b_min; c4_min=c_min;                
        a5_min=a_min; b5_min=b_min; c5_min=c_min;               
        a6_min=a_min; b6_min=b_min; c6_min=c_min;                
        a7_min=a_min; b7_min=b_min; c7_min=c_min;                
        a1_max=inf; b1_max=inf; c1_max=inf;                
        a2_max=inf; b2_max=inf; c2_max=inf;                
        a3_max=inf; b3_max=inf; c3_max=inf;          
        a4_max=inf; b4_max=inf; c4_max=inf;         
        a5_max=inf; b5_max=inf; c5_max=inf;         
        a6_max=inf; b6_max=inf; c6_max=inf;         
        a7_max=inf; b7_max=inf; c7_max=inf; 
        opts.Lower = [a1_min b1_min c1_min a2_min b2_min c2_min a3_min b3_min c3_min a4_min b4_min c4_min a5_min b5_min c5_min a6_min b6_min c6_min a7_min b7_min c7_min];
        opts.Upper = [a1_max b1_max c1_max a2_max b2_max c2_max a3_max b3_max c3_max a4_max b4_max c4_max a5_max b5_max c5_max a6_max b6_max c6_max a7_max b7_max c7_max];
        opts.Robust = 'Bisquare';
        [fitresult7,gof7] = fit(logX,Y_interp,ft,opts);
        y_G7 = fitresult7(logX);
        A1_opt_G7=fitresult7.a1;
        A2_opt_G7=fitresult7.a2;
        A3_opt_G7=fitresult7.a3;
        A4_opt_G7=fitresult7.a4;
        A5_opt_G7=fitresult7.a5;
        A6_opt_G7=fitresult7.a6;
        A7_opt_G7=fitresult7.a7;
        B1_opt_G7=fitresult7.b1;
        B2_opt_G7=fitresult7.b2;
        B3_opt_G7=fitresult7.b3;
        B4_opt_G7=fitresult7.b4;
        B5_opt_G7=fitresult7.b5;
        B6_opt_G7=fitresult7.b6;
        B7_opt_G7=fitresult7.b7;
        C1_opt_G7=fitresult7.c1;
        C2_opt_G7=fitresult7.c2;
        C3_opt_G7=fitresult7.c3;    
        C4_opt_G7=fitresult7.c4;    
        C5_opt_G7=fitresult7.c5;    
        C6_opt_G7=fitresult7.c6;      
        C7_opt_G7=fitresult7.c7;              
        df_a1=exp(-((logX-B1_opt_G7).^2./C1_opt_G7^2));
        df_a2=exp(-((logX-B2_opt_G7).^2./C2_opt_G7^2));
        df_a3=exp(-((logX-B3_opt_G7).^2./C3_opt_G7^2));
        df_a4=exp(-((logX-B4_opt_G7).^2./C4_opt_G7^2));
        df_a5=exp(-((logX-B5_opt_G7).^2./C5_opt_G7^2));
        df_a6=exp(-((logX-B6_opt_G7).^2./C6_opt_G7^2));
        df_a7=exp(-((logX-B7_opt_G7).^2./C7_opt_G7^2));
        df_b1=exp(-((logX-B1_opt_G7).^2./C1_opt_G7^2)) .* 2.*A1_opt_G7.*(logX-B1_opt_G7)./C1_opt_G7^2;
        df_b2=exp(-((logX-B2_opt_G7).^2./C2_opt_G7^2)) .* 2.*A2_opt_G7.*(logX-B2_opt_G7)./C2_opt_G7^2;
        df_b3=exp(-((logX-B3_opt_G7).^2./C3_opt_G7^2)) .* 2.*A3_opt_G7.*(logX-B3_opt_G7)./C3_opt_G7^2;
        df_b4=exp(-((logX-B4_opt_G7).^2./C4_opt_G7^2)) .* 2.*A4_opt_G7.*(logX-B4_opt_G7)./C4_opt_G7^2;
        df_b5=exp(-((logX-B5_opt_G7).^2./C5_opt_G7^2)) .* 2.*A5_opt_G7.*(logX-B5_opt_G7)./C5_opt_G7^2;
        df_b6=exp(-((logX-B6_opt_G7).^2./C6_opt_G7^2)) .* 2.*A6_opt_G7.*(logX-B6_opt_G7)./C6_opt_G7^2;
        df_b7=exp(-((logX-B7_opt_G7).^2./C7_opt_G7^2)) .* 2.*A7_opt_G7.*(logX-B7_opt_G7)./C7_opt_G7^2;
        df_c1=exp(-((logX-B1_opt_G7).^2./C1_opt_G7^2)) .* 2.*A1_opt_G7.*(logX-B1_opt_G7).^2./C1_opt_G7^3;                      
        df_c2=exp(-((logX-B2_opt_G7).^2./C2_opt_G7^2)) .* 2.*A2_opt_G7.*(logX-B2_opt_G7).^2./C2_opt_G7^3;                      
        df_c3=exp(-((logX-B3_opt_G7).^2./C3_opt_G7^2)) .* 2.*A3_opt_G7.*(logX-B3_opt_G7).^2./C3_opt_G7^3;                      
        df_c4=exp(-((logX-B4_opt_G7).^2./C4_opt_G7^2)) .* 2.*A4_opt_G7.*(logX-B4_opt_G7).^2./C4_opt_G7^3;                      
        df_c5=exp(-((logX-B5_opt_G7).^2./C5_opt_G7^2)) .* 2.*A5_opt_G7.*(logX-B5_opt_G7).^2./C5_opt_G7^3;                      
        df_c6=exp(-((logX-B6_opt_G7).^2./C6_opt_G7^2)) .* 2.*A6_opt_G7.*(logX-B6_opt_G7).^2./C6_opt_G7^3;                      
        df_c7=exp(-((logX-B7_opt_G7).^2./C7_opt_G7^2)) .* 2.*A7_opt_G7.*(logX-B7_opt_G7).^2./C7_opt_G7^3;                      
        ci = confint(fitresult7,0.95); ci(isnan(ci))=0;
        SE_a1=(A1_opt_G7-ci(1,1))/1.96;
        SE_b1=(B1_opt_G7-ci(1,2))/1.96;
        SE_c1=(C1_opt_G7-ci(1,3))/1.96;
        SE_a2=(A2_opt_G7-ci(1,4))/1.96;
        SE_b2=(B2_opt_G7-ci(1,5))/1.96;
        SE_c2=(C2_opt_G7-ci(1,6))/1.96;
        SE_a3=(A3_opt_G7-ci(1,7))/1.96;
        SE_b3=(B3_opt_G7-ci(1,8))/1.96;
        SE_c3=(C3_opt_G7-ci(1,9))/1.96;
        SE_a4=(A4_opt_G7-ci(1,10))/1.96;
        SE_b4=(B4_opt_G7-ci(1,11))/1.96;
        SE_c4=(C4_opt_G7-ci(1,12))/1.96;
        SE_a5=(A5_opt_G7-ci(1,13))/1.96;
        SE_b5=(B5_opt_G7-ci(1,14))/1.96;
        SE_c5=(C5_opt_G7-ci(1,15))/1.96;
        SE_a6=(A6_opt_G7-ci(1,16))/1.96;
        SE_b6=(B6_opt_G7-ci(1,17))/1.96;
        SE_c6=(C6_opt_G7-ci(1,18))/1.96;        
        SE_a7=(A7_opt_G7-ci(1,19))/1.96;
        SE_b7=(B7_opt_G7-ci(1,20))/1.96;
        SE_c7=(C7_opt_G7-ci(1,21))/1.96; 
        SE_G7=abs(sqrt( (df_a1.*SE_a1).^2+(df_b1.*SE_b1).^2+(df_c1.*SE_c1).^2 + (df_a2.*SE_a2).^2+(df_b2.*SE_b2).^2+(df_c2.*SE_c2).^2 + (df_a3.*SE_a3).^2+(df_b3.*SE_b3).^2+(df_c3.*SE_c3).^2  + (df_a4.*SE_a4).^2+(df_b4.*SE_b4).^2+(df_c4.*SE_c4).^2  + (df_a5.*SE_a5).^2+(df_b5.*SE_b5).^2+(df_c5.*SE_c5).^2  + (df_a6.*SE_a6).^2+(df_b6.*SE_b6).^2+(df_c6.*SE_c6).^2  + (df_a7.*SE_a7).^2+(df_b7.*SE_b7).^2+(df_c7.*SE_c7).^2    ));
        y_G7_U=y_G7+SE_G7;
        y_G7_L=y_G7-SE_G7;  
        V_G7 = trapz(logX,y_G7); str_G7={['7-modal fit: ',str_Y,'=',num2str(V_G7,str_precision)]}; 
        y_G7_U=y_G7+SE_G7; y_G7_L=y_G7-SE_G7;  
        GMM_coef(7,:)=[7,A1_opt_G7,A2_opt_G7,A3_opt_G7,A4_opt_G7,A5_opt_G7,A6_opt_G7,A7_opt_G7,zeros(1,1),B1_opt_G7,B2_opt_G7,B3_opt_G7,B4_opt_G7,B5_opt_G7,B6_opt_G7,B7_opt_G7,zeros(1,1),C1_opt_G7,C2_opt_G7,C3_opt_G7,C4_opt_G7,C5_opt_G7,C6_opt_G7,C7_opt_G7,zeros(1,1)];
        GMM_se(7,:)=[7,SE_a1,SE_a2,SE_a3,SE_a4,SE_a5,SE_a6,SE_a7,zeros(1,1),SE_b1,SE_b2,SE_b3,SE_b4,SE_b5,SE_b6,SE_b7,zeros(1,1),SE_c1,SE_c2,SE_c3,SE_c4,SE_c5,SE_c6,SE_c7,zeros(1,1)];
        
        if N_interp>24

            % GAUSS-8 FIT
            ft = fittype( 'gauss8' );
            opts = fitoptions( ft ); 
            a1_min=a_min; b1_min=b_min; c1_min=c_min;                    
            a2_min=a_min; b2_min=b_min; c2_min=c_min;                    
            a3_min=a_min; b3_min=b_min; c3_min=c_min;                    
            a4_min=a_min; b4_min=b_min; c4_min=c_min;                    
            a5_min=a_min; b5_min=b_min; c5_min=c_min;                    
            a6_min=a_min; b6_min=b_min; c6_min=c_min;                    
            a7_min=a_min; b7_min=b_min; c7_min=c_min;                    
            a8_min=a_min; b8_min=b_min; c8_min=c_min;                    
            a1_max=inf; b1_max=inf; c1_max=inf;                    
            a2_max=inf; b2_max=inf; c2_max=inf;                    
            a3_max=inf; b3_max=inf; c3_max=inf;              
            a4_max=inf; b4_max=inf; c4_max=inf;             
            a5_max=inf; b5_max=inf; c5_max=inf;             
            a6_max=inf; b6_max=inf; c6_max=inf;             
            a7_max=inf; b7_max=inf; c7_max=inf;             
            a8_max=inf; b8_max=inf; c8_max=inf;             
            opts.Lower = [a1_min b1_min c1_min a2_min b2_min c2_min a3_min b3_min c3_min a4_min b4_min c4_min a5_min b5_min c5_min a6_min b6_min c6_min a7_min b7_min c7_min a8_min b8_min c8_min];
            opts.Upper = [a1_max b1_max c1_max a2_max b2_max c2_max a3_max b3_max c3_max a4_max b4_max c4_max a5_max b5_max c5_max a6_max b6_max c6_max a7_max b7_max c7_max a8_max b8_max c8_max];
            opts.Robust = 'Bisquare';
            [fitresult8,gof8] = fit(logX,Y_interp,ft,opts);
            y_G8= fitresult8(logX);
            A1_opt_G8=fitresult8.a1;
            A2_opt_G8=fitresult8.a2;
            A3_opt_G8=fitresult8.a3;
            A4_opt_G8=fitresult8.a4;
            A5_opt_G8=fitresult8.a5;
            A6_opt_G8=fitresult8.a6;
            A7_opt_G8=fitresult8.a7;
            A8_opt_G8=fitresult8.a8;
            B1_opt_G8=fitresult8.b1;
            B2_opt_G8=fitresult8.b2;
            B3_opt_G8=fitresult8.b3;
            B4_opt_G8=fitresult8.b4;
            B5_opt_G8=fitresult8.b5;
            B6_opt_G8=fitresult8.b6;
            B7_opt_G8=fitresult8.b7;
            B8_opt_G8=fitresult8.b8;
            C1_opt_G8=fitresult8.c1;
            C2_opt_G8=fitresult8.c2;
            C3_opt_G8=fitresult8.c3;    
            C4_opt_G8=fitresult8.c4;    
            C5_opt_G8=fitresult8.c5;    
            C6_opt_G8=fitresult8.c6;      
            C7_opt_G8=fitresult8.c7;      
            C8_opt_G8=fitresult8.c8;                      
            df_a1=exp(-((logX-B1_opt_G8).^2./C1_opt_G8^2));
            df_a2=exp(-((logX-B2_opt_G8).^2./C2_opt_G8^2));
            df_a3=exp(-((logX-B3_opt_G8).^2./C3_opt_G8^2));
            df_a4=exp(-((logX-B4_opt_G8).^2./C4_opt_G8^2));
            df_a5=exp(-((logX-B5_opt_G8).^2./C5_opt_G8^2));
            df_a6=exp(-((logX-B6_opt_G8).^2./C6_opt_G8^2));
            df_a7=exp(-((logX-B7_opt_G8).^2./C7_opt_G8^2));
            df_a8=exp(-((logX-B8_opt_G8).^2./C8_opt_G8^2));
            df_b1=exp(-((logX-B1_opt_G8).^2./C1_opt_G8^2)) .* 2.*A1_opt_G8.*(logX-B1_opt_G8)./C1_opt_G8^2;
            df_b2=exp(-((logX-B2_opt_G8).^2./C2_opt_G8^2)) .* 2.*A2_opt_G8.*(logX-B2_opt_G8)./C2_opt_G8^2;
            df_b3=exp(-((logX-B3_opt_G8).^2./C3_opt_G8^2)) .* 2.*A3_opt_G8.*(logX-B3_opt_G8)./C3_opt_G8^2;
            df_b4=exp(-((logX-B4_opt_G8).^2./C4_opt_G8^2)) .* 2.*A4_opt_G8.*(logX-B4_opt_G8)./C4_opt_G8^2;
            df_b5=exp(-((logX-B5_opt_G8).^2./C5_opt_G8^2)) .* 2.*A5_opt_G8.*(logX-B5_opt_G8)./C5_opt_G8^2;
            df_b6=exp(-((logX-B6_opt_G8).^2./C6_opt_G8^2)) .* 2.*A6_opt_G8.*(logX-B6_opt_G8)./C6_opt_G8^2;
            df_b7=exp(-((logX-B7_opt_G8).^2./C7_opt_G8^2)) .* 2.*A7_opt_G8.*(logX-B7_opt_G8)./C7_opt_G8^2;
            df_b8=exp(-((logX-B8_opt_G8).^2./C8_opt_G8^2)) .* 2.*A8_opt_G8.*(logX-B8_opt_G8)./C8_opt_G8^2;
            df_c1=exp(-((logX-B1_opt_G8).^2./C1_opt_G8^2)) .* 2.*A1_opt_G8.*(logX-B1_opt_G8).^2./C1_opt_G8^3;                      
            df_c2=exp(-((logX-B2_opt_G8).^2./C2_opt_G8^2)) .* 2.*A2_opt_G8.*(logX-B2_opt_G8).^2./C2_opt_G8^3;                      
            df_c3=exp(-((logX-B3_opt_G8).^2./C3_opt_G8^2)) .* 2.*A3_opt_G8.*(logX-B3_opt_G8).^2./C3_opt_G8^3;                      
            df_c4=exp(-((logX-B4_opt_G8).^2./C4_opt_G8^2)) .* 2.*A4_opt_G8.*(logX-B4_opt_G8).^2./C4_opt_G8^3;                      
            df_c5=exp(-((logX-B5_opt_G8).^2./C5_opt_G8^2)) .* 2.*A5_opt_G8.*(logX-B5_opt_G8).^2./C5_opt_G8^3;                      
            df_c6=exp(-((logX-B6_opt_G8).^2./C6_opt_G8^2)) .* 2.*A6_opt_G8.*(logX-B6_opt_G8).^2./C6_opt_G8^3;                      
            df_c7=exp(-((logX-B7_opt_G8).^2./C7_opt_G8^2)) .* 2.*A7_opt_G8.*(logX-B7_opt_G8).^2./C7_opt_G8^3;                      
            df_c8=exp(-((logX-B8_opt_G8).^2./C8_opt_G8^2)) .* 2.*A8_opt_G8.*(logX-B8_opt_G8).^2./C8_opt_G8^3;                      
            ci = confint(fitresult8,0.95); ci(isnan(ci))=0;
            SE_a1=(A1_opt_G8-ci(1,1))/1.96;
            SE_b1=(B1_opt_G8-ci(1,2))/1.96;
            SE_c1=(C1_opt_G8-ci(1,3))/1.96;
            SE_a2=(A2_opt_G8-ci(1,4))/1.96;
            SE_b2=(B2_opt_G8-ci(1,5))/1.96;
            SE_c2=(C2_opt_G8-ci(1,6))/1.96;
            SE_a3=(A3_opt_G8-ci(1,7))/1.96;
            SE_b3=(B3_opt_G8-ci(1,8))/1.96;
            SE_c3=(C3_opt_G8-ci(1,9))/1.96;
            SE_a4=(A4_opt_G8-ci(1,10))/1.96;
            SE_b4=(B4_opt_G8-ci(1,11))/1.96;
            SE_c4=(C4_opt_G8-ci(1,12))/1.96;
            SE_a5=(A5_opt_G8-ci(1,13))/1.96;
            SE_b5=(B5_opt_G8-ci(1,14))/1.96;
            SE_c5=(C5_opt_G8-ci(1,15))/1.96;
            SE_a6=(A6_opt_G8-ci(1,16))/1.96;
            SE_b6=(B6_opt_G8-ci(1,17))/1.96;
            SE_c6=(C6_opt_G8-ci(1,18))/1.96;        
            SE_a7=(A7_opt_G8-ci(1,19))/1.96;
            SE_b7=(B7_opt_G8-ci(1,20))/1.96;
            SE_c7=(C7_opt_G8-ci(1,21))/1.96; 
            SE_a8=(A8_opt_G8-ci(1,22))/1.96;
            SE_b8=(B8_opt_G8-ci(1,23))/1.96;
            SE_c8=(C8_opt_G8-ci(1,24))/1.96;         
            SE_G8=abs(sqrt( (df_a1.*SE_a1).^2+(df_b1.*SE_b1).^2+(df_c1.*SE_c1).^2 + (df_a2.*SE_a2).^2+(df_b2.*SE_b2).^2+(df_c2.*SE_c2).^2 + (df_a3.*SE_a3).^2+(df_b3.*SE_b3).^2+(df_c3.*SE_c3).^2  + (df_a4.*SE_a4).^2+(df_b4.*SE_b4).^2+(df_c4.*SE_c4).^2  + (df_a5.*SE_a5).^2+(df_b5.*SE_b5).^2+(df_c5.*SE_c5).^2  + (df_a6.*SE_a6).^2+(df_b6.*SE_b6).^2+(df_c6.*SE_c6).^2  + (df_a7.*SE_a7).^2+(df_b7.*SE_b7).^2+(df_c7.*SE_c7).^2  + (df_a8.*SE_a8).^2+(df_b8.*SE_b8).^2+(df_c8.*SE_c8).^2    ));
            y_G8_U=y_G8+SE_G8;
            y_G8_L=y_G8-SE_G8;  
            V_G8 = trapz(logX,y_G8); str_G8={['8-modal fit: ',str_Y,'=',num2str(V_G8,str_precision)]};
            y_G8_U=y_G8+SE_G8; y_G8_L=y_G8-SE_G8; 
            GMM_coef(8,:)=[8,A1_opt_G8,A2_opt_G8,A3_opt_G8,A4_opt_G8,A5_opt_G8,A6_opt_G8,A7_opt_G8,A8_opt_G8,B1_opt_G8,B2_opt_G8,B3_opt_G8,B4_opt_G8,B5_opt_G8,B6_opt_G8,B7_opt_G8,B8_opt_G8,C1_opt_G8,C2_opt_G8,C3_opt_G8,C4_opt_G8,C5_opt_G8,C6_opt_G8,C7_opt_G8,C8_opt_G8];
            GMM_se(8,:)=[8,SE_a1,SE_a2,SE_a3,SE_a4,SE_a5,SE_a6,SE_a7,SE_a8,SE_b1,SE_b2,SE_b3,SE_b4,SE_b5,SE_b6,SE_b7,SE_b8,SE_c1,SE_c2,SE_c3,SE_c4,SE_c5,SE_c6,SE_c7,SE_c8];
        else
            GMM_coef(8,:)=[8,zeros(1,24)];
            GMM_se(8,:)=[8,zeros(1,24)];
        end

        if isequal(FLAG_plot_GMM,0)
            % n=1
            A1=fitresult1.a1; % V
            B1=fitresult1.b1; % ln(r)
            C1=fitresult1.c1; % ln(sd)
            G1_1    = A1*exp(-((logX-B1)/C1).^2); % uni-modal log-normal distribution
            V_G1_1  = trapz(logX,G1_1);
            r_G1_1  = exp(B1);
            sd_G1_1 = (1/sqrt(2))*(C1); 
            % n=2            
            A1=fitresult2.a1;
            A2=fitresult2.a2;
            B1=fitresult2.b1;
            B2=fitresult2.b2;
            C1=fitresult2.c1;
            C2=fitresult2.c2;
            G2_1    = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2); % bi-modal log-normal distribution
            G2_2    = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2); % bi-modal log-normal distribution
            V_G2_1  = trapz(logX,G2_1);
            V_G2_2  = trapz(logX,G2_2);
            r_G2_1  = exp(B1);
            r_G2_2  = exp(B2);
            sd_G2_1 = (1/sqrt(2))*(C1);    
            sd_G2_2 = (1/sqrt(2))*(C2); 
            % n=3            
            A1=fitresult3.a1;
            A2=fitresult3.a2;
            A3=fitresult3.a3;
            B1=fitresult3.b1;
            B2=fitresult3.b2;
            B3=fitresult3.b3;
            C1=abs(fitresult3.c1);
            C2=abs(fitresult3.c2);
            C3=abs(fitresult3.c3);
            G3_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2); % tri-modal log-normal distribution
            G3_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2); % tri-modal log-normal distribution
            G3_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2); % tri-modal log-normal distribution
            V_G3_1 = trapz(logX,G3_1);
            V_G3_2 = trapz(logX,G3_2);
            V_G3_3 = trapz(logX,G3_3);
            r_G3_1 = exp(B1);
            r_G3_2 = exp(B2);
            r_G3_3 = exp(B3);
            sd_G3_1 = (1/sqrt(2))*(C1);    
            sd_G3_2 = (1/sqrt(2))*(C2);    
            sd_G3_3 = (1/sqrt(2))*(C3); 
            % n=4            
            A1=fitresult4.a1;
            A2=fitresult4.a2;
            A3=fitresult4.a3;
            A4=fitresult4.a4;
            B1=fitresult4.b1;
            B2=fitresult4.b2;
            B3=fitresult4.b3;
            B4=fitresult4.b4;
            C1=fitresult4.c1;
            C2=fitresult4.c2;
            C3=fitresult4.c3;
            C4=fitresult4.c4;
            G4_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2); % quad-modal log-normal distribution
            G4_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2); % quad-modal log-normal distribution
            G4_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2); % quad-modal log-normal distribution
            G4_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2); % quad-modal log-normal distribution
            V_G4_1 = trapz(logX,G4_1);
            V_G4_2 = trapz(logX,G4_2);
            V_G4_3 = trapz(logX,G4_3);
            V_G4_4 = trapz(logX,G4_4);
            r_G4_1 = exp(B1);
            r_G4_2 = exp(B2);
            r_G4_3 = exp(B3);
            r_G4_4 = exp(B4);
            sd_G4_1 = (1/sqrt(2))*(C1);    
            sd_G4_2 = (1/sqrt(2))*(C2);    
            sd_G4_3 = (1/sqrt(2))*(C3);    
            sd_G4_4 = (1/sqrt(2))*(C4);  
            % n=5            
            A1=fitresult5.a1;
            A2=fitresult5.a2;
            A3=fitresult5.a3;
            A4=fitresult5.a4;
            A5=fitresult5.a5;
            B1=fitresult5.b1;
            B2=fitresult5.b2;
            B3=fitresult5.b3;
            B4=fitresult5.b4;
            B5=fitresult5.b5;
            C1=fitresult5.c1;
            C2=fitresult5.c2;
            C3=fitresult5.c3;
            C4=fitresult5.c4;
            C5=fitresult5.c5;
            G5_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            G5_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            G5_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            G5_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            G5_5 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+A5*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            V_G5_1 = trapz(logX,G5_1);
            V_G5_2 = trapz(logX,G5_2);
            V_G5_3 = trapz(logX,G5_3);
            V_G5_4 = trapz(logX,G5_4);
            V_G5_5 = trapz(logX,G5_5);
            r_G5_1 = exp(B1);
            r_G5_2 = exp(B2);
            r_G5_3 = exp(B3);
            r_G5_4 = exp(B4);
            r_G5_5 = exp(B5);
            sd_G5_1 = (1/sqrt(2))*(C1);    
            sd_G5_2 = (1/sqrt(2))*(C2);    
            sd_G5_3 = (1/sqrt(2))*(C3);    
            sd_G5_4 = (1/sqrt(2))*(C4);    
            sd_G5_5 = (1/sqrt(2))*(C5); 
            % n=6            
            A1=fitresult6.a1;
            A2=fitresult6.a2;
            A3=fitresult6.a3;
            A4=fitresult6.a4;
            A5=fitresult6.a5;
            A6=fitresult6.a6;
            B1=fitresult6.b1;
            B2=fitresult6.b2;
            B3=fitresult6.b3;
            B4=fitresult6.b4;
            B5=fitresult6.b5;
            B6=fitresult6.b6;
            C1=fitresult6.c1;
            C2=fitresult6.c2;
            C3=fitresult6.c3;
            C4=fitresult6.c4;
            C5=fitresult6.c5;
            C6=fitresult6.c6;
            G6_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_5 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+A5*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_6 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+A6*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            V_G6_1 = trapz(logX,G6_1);
            V_G6_2 = trapz(logX,G6_2);
            V_G6_3 = trapz(logX,G6_3);
            V_G6_4 = trapz(logX,G6_4);
            V_G6_5 = trapz(logX,G6_5);
            V_G6_6 = trapz(logX,G6_6);
            r_G6_1 = exp(B1);
            r_G6_2 = exp(B2);
            r_G6_3 = exp(B3);
            r_G6_4 = exp(B4);
            r_G6_5 = exp(B5);
            r_G6_6 = exp(B6);
            sd_G6_1 = (1/sqrt(2))*(C1);    
            sd_G6_2 = (1/sqrt(2))*(C2);    
            sd_G6_3 = (1/sqrt(2))*(C3);    
            sd_G6_4 = (1/sqrt(2))*(C4);    
            sd_G6_5 = (1/sqrt(2))*(C5);    
            sd_G6_6 = (1/sqrt(2))*(C6);  
            % n=7            
            A1=fitresult7.a1;
            A2=fitresult7.a2;
            A3=fitresult7.a3;
            A4=fitresult7.a4;
            A5=fitresult7.a5;
            A6=fitresult7.a6;
            A7=fitresult7.a7;
            B1=fitresult7.b1;
            B2=fitresult7.b2;
            B3=fitresult7.b3;
            B4=fitresult7.b4;
            B5=fitresult7.b5;
            B6=fitresult7.b6;
            B7=fitresult7.b7;
            C1=fitresult7.c1;
            C2=fitresult7.c2;
            C3=fitresult7.c3;
            C4=fitresult7.c4;
            C5=fitresult7.c5;
            C6=fitresult7.c6;
            C7=fitresult7.c7;
            G7_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_5 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+A5*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_6 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+A6*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_7 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+A7*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            V_G7_1 = trapz(logX,G7_1);
            V_G7_2 = trapz(logX,G7_2);
            V_G7_3 = trapz(logX,G7_3);
            V_G7_4 = trapz(logX,G7_4);
            V_G7_5 = trapz(logX,G7_5);
            V_G7_6 = trapz(logX,G7_6);
            V_G7_7 = trapz(logX,G7_7);
            r_G7_1 = exp(B1);
            r_G7_2 = exp(B2);
            r_G7_3 = exp(B3);
            r_G7_4 = exp(B4);
            r_G7_5 = exp(B5);
            r_G7_6 = exp(B6);
            r_G7_7 = exp(B7);
            sd_G7_1 = (1/sqrt(2))*(C1);    
            sd_G7_2 = (1/sqrt(2))*(C2);    
            sd_G7_3 = (1/sqrt(2))*(C3);    
            sd_G7_4 = (1/sqrt(2))*(C4);    
            sd_G7_5 = (1/sqrt(2))*(C5);    
            sd_G7_6 = (1/sqrt(2))*(C6);    
            sd_G7_7 = (1/sqrt(2))*(C7); 
            % n=8            
            A1=fitresult8.a1;
            A2=fitresult8.a2;
            A3=fitresult8.a3;
            A4=fitresult8.a4;
            A5=fitresult8.a5;
            A6=fitresult8.a6;
            A7=fitresult8.a7;
            A8=fitresult8.a8;
            B1=fitresult8.b1;
            B2=fitresult8.b2;
            B3=fitresult8.b3;
            B4=fitresult8.b4;
            B5=fitresult8.b5;
            B6=fitresult8.b6;
            B7=fitresult8.b7;
            B8=fitresult8.b8;
            C1=fitresult8.c1;
            C2=fitresult8.c2;
            C3=fitresult8.c3;
            C4=fitresult8.c4;
            C5=fitresult8.c5;
            C6=fitresult8.c6;
            C7=fitresult8.c7;
            C8=fitresult8.c8;
            G8_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
            G8_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
            G8_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
            G8_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
            G8_5 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+A5*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
            G8_6 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+A6*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
            G8_7 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+A7*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
            G8_8 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+A8*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
            V_G8_1 = trapz(logX,G8_1);
            V_G8_2 = trapz(logX,G8_2);
            V_G8_3 = trapz(logX,G8_3);
            V_G8_4 = trapz(logX,G8_4);
            V_G8_5 = trapz(logX,G8_5);
            V_G8_6 = trapz(logX,G8_6);
            V_G8_7 = trapz(logX,G8_7);
            V_G8_8 = trapz(logX,G8_8);
            r_G8_1 = exp(B1);
            r_G8_2 = exp(B2);
            r_G8_3 = exp(B3);
            r_G8_4 = exp(B4);
            r_G8_5 = exp(B5);
            r_G8_6 = exp(B6);
            r_G8_7 = exp(B7);
            r_G8_8 = exp(B8);
            sd_G8_1 = (1/sqrt(2))*(C1);    
            sd_G8_2 = (1/sqrt(2))*(C2);    
            sd_G8_3 = (1/sqrt(2))*(C3);    
            sd_G8_4 = (1/sqrt(2))*(C4);    
            sd_G8_5 = (1/sqrt(2))*(C5);    
            sd_G8_6 = (1/sqrt(2))*(C6);    
            sd_G8_7 = (1/sqrt(2))*(C7);    
            sd_G8_8 = (1/sqrt(2))*(C8);           
        else
            str_AERONET={['Raw data: ',str_Y,'=',num2str(V_trapz,str_precision)]};        
            if isequal(FLAG_plot_visible,0)
                figure; set(gcf, 'color','white', 'visible','off','units','normalized','outerposition',[0 0 1 1]);
            else   
                figure; set(gcf, 'color','white', 'visible','on','units','normalized','outerposition',[0 0 1 1]);    
            end
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
            hold on
            h_G1 = semilogx(X_interp,y_G1,'color',bright,'LineWidth',2*line_width);
            hold on
            A1=fitresult1.a1; % V
            B1=fitresult1.b1; % ln(r)
            C1=fitresult1.c1; % ln(sd)
            G1_1    = A1*exp(-((logX-B1)/C1).^2); % uni-modal log-normal distribution
            V_G1_1  = trapz(logX,G1_1);
            r_G1_1  = exp(B1);
            sd_G1_1 = (1/sqrt(2))*(C1);  
            x=Y_interp; y=y_G1; n_coeff=3;        
            [res_1,sse_1,R2_1,rmse_1]=my_stats(x,y); b_1=sum(x-y)/n_interp; df_1=n_interp-n_coeff; R2adj_1=1-(1-R2_1)*(n_interp-1)/(df_1-1);                                
            str_1={[str_Y,'(1)=',num2str(V_G1_1,str_precision),' r(1)=',num2str(r_G1_1,str_precision),' \sigma(1)=',num2str(sd_G1_1,str_precision)]};
            h_1=semilogx(X_interp,G1_1,'color',bright,'LineWidth',line_width,'LineStyle','-');
            hold on
            semilogx(X_interp,y_G1_U,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on
            semilogx(X_interp,y_G1_L,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on
            line([exp(B1) exp(B1)],[0 max(G1_1)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B1-sqrt(2*log(2))*sd_G1_1) exp(B1+sqrt(2*log(2))*sd_G1_1)],[max(G1_1)/2 max(G1_1)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on    
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
%             xmin=min(X); xmax=max(X);ymin=0; ymax=max([max(Y_interp),max(Y),max(y_G1)]); axis([xmin xmax ymin ymax]);
            xlabel(str_x,'FontSize',font_size,'FontWeight','bold');ylabel(str_y,'FontSize',font_size,'FontWeight','bold');
            title(['1-modal: ','b=',num2str(b_1,str_precision),' s=',num2str(rmse_1,str_precision),' R^2=',num2str(R2adj_1,str_precision)],'FontSize',font_size,'FontWeight','bold');
            legendText=[];  
            legendText=[legendText,str_G1,str_1];
            legend([h_G1 h_1],legendText(:),'Location','EastOutside','FontSize',font_size);
            ylimits = get(gca,'ylim'); set(gca,'ylim',[0,ylimits(2)]);
            set(gcf,'PaperPositionMode','auto');  
            set(gca,'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XColor','k','YColor','k','LineWidth',line_width,'FontSize',font_size);
            if isequal(FLAG_plot_save,1)
                file_name = fullfile(path_save,['GMM1']); 
                print('-djpeg','-r200',file_name);                                                       
                close              
            end
            
            if isequal(FLAG_plot_visible,0)
                figure; set(gcf, 'color','white', 'visible','off','units','normalized','outerposition',[0 0 1 1]);
            else   
                figure; set(gcf, 'color','white', 'visible','on','units','normalized','outerposition',[0 0 1 1]);    
            end            
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
            hold on     
            h_G2 = semilogx(X_interp,y_G2,'color',bright,'LineWidth',2*line_width);
            hold on
            A1=fitresult2.a1;
            A2=fitresult2.a2;
            B1=fitresult2.b1;
            B2=fitresult2.b2;
            C1=fitresult2.c1;
            C2=fitresult2.c2;
            G2_1    = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2); % bi-modal log-normal distribution
            G2_2    = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2); % bi-modal log-normal distribution
            V_G2_1  = trapz(logX,G2_1);
            V_G2_2  = trapz(logX,G2_2);
            r_G2_1  = exp(B1);
            r_G2_2  = exp(B2);
            sd_G2_1 = (1/sqrt(2))*(C1);    
            sd_G2_2 = (1/sqrt(2))*(C2);    
            x=Y_interp; y=y_G2; n_coeff=6;        
            [res_2,sse_2,R2_2,rmse_2]=my_stats(x,y); b_2=sum(x-y)/n_interp; df_2=n_interp-n_coeff; R2adj_2=1-(1-R2_2)*(n_interp-1)/(df_2-1);                                
            str_1={[str_Y,'(1)=',num2str(V_G2_1,str_precision),' r(1)=',num2str(r_G2_1,str_precision),' \sigma(1)=',num2str(sd_G2_1,str_precision)]};
            str_2={[str_Y,'(2)=',num2str(V_G2_2,str_precision),' r(2)=',num2str(r_G2_2,str_precision),' \sigma(2)=',num2str(sd_G2_2,str_precision)]};
            h_1=semilogx(X_interp,G2_1,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_2=semilogx(X_interp,G2_2,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            semilogx(X_interp,y_G2_U,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on
            semilogx(X_interp,y_G2_L,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on        
            line([exp(B1) exp(B1)],[0 max(G2_1)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2) exp(B2)],[0 max(G2_2)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B1-sqrt(2*log(2))*sd_G2_1) exp(B1+sqrt(2*log(2))*sd_G2_1)],[max(G2_1)/2 max(G2_1)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2-sqrt(2*log(2))*sd_G2_2) exp(B2+sqrt(2*log(2))*sd_G2_2)],[max(G2_2)/2 max(G2_2)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
            xlabel(str_x,'FontSize',font_size,'FontWeight','bold');ylabel(str_y,'FontSize',font_size,'FontWeight','bold');
            title(['2-modal: ','b=',num2str(b_2,str_precision),' s=',num2str(rmse_2,str_precision),' R^2=',num2str(R2adj_2,str_precision)],'FontSize',font_size,'FontWeight','bold');
            legendText=[];  
            legendText=[legendText,str_G2,str_1,str_2];
            legend([h_G2,h_1,h_2],legendText(:),'Location','EastOutside','FontSize',font_size);            
            ylimits = get(gca,'ylim'); set(gca,'ylim',[0,ylimits(2)]);
            set(gcf,'PaperPositionMode','auto');  
            set(gca,'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XColor','k','YColor','k','LineWidth',line_width,'FontSize',font_size);
            if isequal(FLAG_plot_save,1)
                file_name = fullfile(path_save,['GMM2']); 
                print('-djpeg','-r200',file_name);                                                       
                close              
            end
                        
            if isequal(FLAG_plot_visible,0)
                figure; set(gcf, 'color','white', 'visible','off','units','normalized','outerposition',[0 0 1 1]);
            else   
                figure; set(gcf, 'color','white', 'visible','on','units','normalized','outerposition',[0 0 1 1]);    
            end        
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
            hold on      
            h_G3 = semilogx(X_interp,y_G3,'color',bright,'LineWidth',2*line_width);
            hold on
            A1=fitresult3.a1;
            A2=fitresult3.a2;
            A3=fitresult3.a3;
            B1=fitresult3.b1;
            B2=fitresult3.b2;
            B3=fitresult3.b3;
            C1=abs(fitresult3.c1);
            C2=abs(fitresult3.c2);
            C3=abs(fitresult3.c3);
            G3_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2); % tri-modal log-normal distribution
            G3_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2); % tri-modal log-normal distribution
            G3_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2); % tri-modal log-normal distribution
            V_G3_1 = trapz(logX,G3_1);
            V_G3_2 = trapz(logX,G3_2);
            V_G3_3 = trapz(logX,G3_3);
            r_G3_1 = exp(B1);
            r_G3_2 = exp(B2);
            r_G3_3 = exp(B3);
            sd_G3_1 = (1/sqrt(2))*(C1);    
            sd_G3_2 = (1/sqrt(2))*(C2);    
            sd_G3_3 = (1/sqrt(2))*(C3);    
            x=Y_interp; y=y_G3; n_coeff=9;        
            [res_3,sse_3,R2_3,rmse_3]=my_stats(x,y); b_3=sum(x-y)/n_interp; df_3=n_interp-n_coeff; R2adj_3=1-(1-R2_3)*(n_interp-1)/(df_3-1);                                
            str_1={[str_Y,'(1)=',num2str(V_G3_1,str_precision),' r(1)=',num2str(r_G3_1,str_precision),' \sigma(1)=',num2str(sd_G3_1,str_precision)]};
            str_2={[str_Y,'(2)=',num2str(V_G3_2,str_precision),' r(2)=',num2str(r_G3_2,str_precision),' \sigma(2)=',num2str(sd_G3_2,str_precision)]};
            str_3={[str_Y,'(3)=',num2str(V_G3_3,str_precision),' r(3)=',num2str(r_G3_3,str_precision),' \sigma(3)=',num2str(sd_G3_3,str_precision)]};
            h_1=semilogx(X_interp,G3_1,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_2=semilogx(X_interp,G3_2,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_3=semilogx(X_interp,G3_3,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            semilogx(X_interp,y_G3_U,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on
            semilogx(X_interp,y_G3_L,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on           
            line([exp(B1) exp(B1)],[0 max(G3_1)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2) exp(B2)],[0 max(G3_2)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3) exp(B3)],[0 max(G3_3)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B1-sqrt(2*log(2))*sd_G3_1) exp(B1+sqrt(2*log(2))*sd_G3_1)],[max(G3_1)/2 max(G3_1)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2-sqrt(2*log(2))*sd_G3_2) exp(B2+sqrt(2*log(2))*sd_G3_2)],[max(G3_2)/2 max(G3_2)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3-sqrt(2*log(2))*sd_G3_3) exp(B3+sqrt(2*log(2))*sd_G3_3)],[max(G3_3)/2 max(G3_3)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey);   
            xlabel(str_x,'FontSize',font_size,'FontWeight','bold');ylabel(str_y,'FontSize',font_size,'FontWeight','bold');
            title(['3-modal: ','b=',num2str(b_3,str_precision),' s=',num2str(rmse_3,str_precision),' R^2=',num2str(R2adj_3,str_precision)],'FontSize',font_size,'FontWeight','bold');
            legendText=[];  
            legendText=[legendText,str_G3,str_1,str_2,str_3];
            legend([h_G3,h_1,h_2,h_3],legendText(:),'Location','EastOutside','FontSize',font_size);
            ylimits = get(gca,'ylim'); set(gca,'ylim',[0,ylimits(2)]);
            set(gcf,'PaperPositionMode','auto');  
            set(gca,'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XColor','k','YColor','k','LineWidth',line_width,'FontSize',font_size);
            if isequal(FLAG_plot_save,1)
                file_name = fullfile(path_save,['GMM3']); 
                print('-djpeg','-r200',file_name);                                                       
                close              
            end
            
            if isequal(FLAG_plot_visible,0)                 
                figure; set(gcf, 'color','white', 'visible','off','units','normalized','outerposition',[0 0 1 1]);             
            else
                figure; set(gcf, 'color','white', 'visible','on','units','normalized','outerposition',[0 0 1 1]);                 
            end;            
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
            hold on    
            h_G4 = semilogx(X_interp,y_G4,'color',bright,'LineWidth',2*line_width);
            hold on
            A1=fitresult4.a1;
            A2=fitresult4.a2;
            A3=fitresult4.a3;
            A4=fitresult4.a4;
            B1=fitresult4.b1;
            B2=fitresult4.b2;
            B3=fitresult4.b3;
            B4=fitresult4.b4;
            C1=fitresult4.c1;
            C2=fitresult4.c2;
            C3=fitresult4.c3;
            C4=fitresult4.c4;
            G4_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2); % quad-modal log-normal distribution
            G4_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2); % quad-modal log-normal distribution
            G4_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2); % quad-modal log-normal distribution
            G4_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2); % quad-modal log-normal distribution
            V_G4_1 = trapz(logX,G4_1);
            V_G4_2 = trapz(logX,G4_2);
            V_G4_3 = trapz(logX,G4_3);
            V_G4_4 = trapz(logX,G4_4);
            r_G4_1 = exp(B1);
            r_G4_2 = exp(B2);
            r_G4_3 = exp(B3);
            r_G4_4 = exp(B4);
            sd_G4_1 = (1/sqrt(2))*(C1);    
            sd_G4_2 = (1/sqrt(2))*(C2);    
            sd_G4_3 = (1/sqrt(2))*(C3);    
            sd_G4_4 = (1/sqrt(2))*(C4);    
            x=Y_interp; y=y_G4;  n_coeff=12;
            [res_4,sse_4,R2_4,rmse_4]=my_stats(x,y); b_4=sum(x-y)/n_interp; df_4=n_interp-n_coeff; R2adj_4=1-(1-R2_4)*(n_interp-1)/(df_4-1);              
            str_1={[str_Y,'(1)=',num2str(V_G4_1,str_precision),' r(1)=',num2str(r_G4_1,str_precision),' \sigma(1)=',num2str(sd_G4_1,str_precision)]};
            str_2={[str_Y,'(2)=',num2str(V_G4_2,str_precision),' r(2)=',num2str(r_G4_2,str_precision),' \sigma(2)=',num2str(sd_G4_2,str_precision)]}; 
            str_3={[str_Y,'(3)=',num2str(V_G4_3,str_precision),' r(3)=',num2str(r_G4_3,str_precision),' \sigma(3)=',num2str(sd_G4_3,str_precision)]};
            str_4={[str_Y,'(4)=',num2str(V_G4_4,str_precision),' r(4)=',num2str(r_G4_4,str_precision),' \sigma(4)=',num2str(sd_G4_4,str_precision)]};
            h_1=semilogx(X_interp,G4_1,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_2=semilogx(X_interp,G4_2,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_3=semilogx(X_interp,G4_3,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_4=semilogx(X_interp,G4_4,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            semilogx(X_interp,y_G4_U,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on
            semilogx(X_interp,y_G4_L,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on           
            line([exp(B1) exp(B1)],[0 max(G4_1)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2) exp(B2)],[0 max(G4_2)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3) exp(B3)],[0 max(G4_3)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B4) exp(B4)],[0 max(G4_4)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B1-sqrt(2*log(2))*sd_G4_1) exp(B1+sqrt(2*log(2))*sd_G4_1)],[max(G4_1)/2 max(G4_1)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2-sqrt(2*log(2))*sd_G4_2) exp(B2+sqrt(2*log(2))*sd_G4_2)],[max(G4_2)/2 max(G4_2)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3-sqrt(2*log(2))*sd_G4_3) exp(B3+sqrt(2*log(2))*sd_G4_3)],[max(G4_3)/2 max(G4_3)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B4-sqrt(2*log(2))*sd_G4_4) exp(B4+sqrt(2*log(2))*sd_G4_4)],[max(G4_4)/2 max(G4_4)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey);   
            xlabel(str_x,'FontSize',font_size,'FontWeight','bold');ylabel(str_y,'FontSize',font_size,'FontWeight','bold');
            title(['4-modal: ','b=',num2str(b_4,str_precision),' s=',num2str(rmse_4,str_precision),' R^2=',num2str(R2adj_4,str_precision)],'FontSize',font_size,'FontWeight','bold');
            legendText=[];  
            legendText=[legendText,str_G4,str_1,str_2,str_3,str_4];
            legend([h_G4,h_1,h_2,h_3,h_4],legendText(:),'Location','EastOutside','FontSize',font_size);
            ylimits = get(gca,'ylim'); set(gca,'ylim',[0,ylimits(2)]);
            set(gcf,'PaperPositionMode','auto');  
            set(gca,'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XColor','k','YColor','k','LineWidth',line_width,'FontSize',font_size);
            if isequal(FLAG_plot_save,1)
                file_name = fullfile(path_save,['GMM4']); 
                print('-djpeg','-r200',file_name);                                                       
                close              
            end

            if isequal(FLAG_plot_visible,0)
                figure; set(gcf, 'color','white', 'visible','off','units','normalized','outerposition',[0 0 1 1]);
            else   
                figure; set(gcf, 'color','white', 'visible','on','units','normalized','outerposition',[0 0 1 1]);    
            end        
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
            hold on  
            h_G5 = semilogx(X_interp,y_G5,'color',bright,'LineWidth',2*line_width);
            hold on
            A1=fitresult5.a1;
            A2=fitresult5.a2;
            A3=fitresult5.a3;
            A4=fitresult5.a4;
            A5=fitresult5.a5;
            B1=fitresult5.b1;
            B2=fitresult5.b2;
            B3=fitresult5.b3;
            B4=fitresult5.b4;
            B5=fitresult5.b5;
            C1=fitresult5.c1;
            C2=fitresult5.c2;
            C3=fitresult5.c3;
            C4=fitresult5.c4;
            C5=fitresult5.c5;
            G5_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            G5_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            G5_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            G5_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            G5_5 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+A5*exp(-((logX-B5)/C5).^2); % 5-modal log-normal distribution
            V_G5_1 = trapz(logX,G5_1);
            V_G5_2 = trapz(logX,G5_2);
            V_G5_3 = trapz(logX,G5_3);
            V_G5_4 = trapz(logX,G5_4);
            V_G5_5 = trapz(logX,G5_5);
            r_G5_1 = exp(B1);
            r_G5_2 = exp(B2);
            r_G5_3 = exp(B3);
            r_G5_4 = exp(B4);
            r_G5_5 = exp(B5);
            sd_G5_1 = (1/sqrt(2))*(C1);    
            sd_G5_2 = (1/sqrt(2))*(C2);    
            sd_G5_3 = (1/sqrt(2))*(C3);    
            sd_G5_4 = (1/sqrt(2))*(C4);    
            sd_G5_5 = (1/sqrt(2))*(C5);    
            x=Y_interp; y=y_G5;  n_coeff=15;
            [res_5,sse_5,R2_5,rmse_5]=my_stats(x,y); b_5=sum(x-y)/n_interp; df_5=n_interp-n_coeff; R2adj_5=1-(1-R2_5)*(n_interp-1)/(df_5-1);              
            str_1={[str_Y,'(1)=',num2str(V_G5_1,str_precision),' r(1)=',num2str(r_G5_1,str_precision),' \sigma(1)=',num2str(sd_G5_1,str_precision)]};
            str_2={[str_Y,'(2)=',num2str(V_G5_2,str_precision),' r(2)=',num2str(r_G5_2,str_precision),' \sigma(2)=',num2str(sd_G5_2,str_precision)]}; 
            str_3={[str_Y,'(3)=',num2str(V_G5_3,str_precision),' r(3)=',num2str(r_G5_3,str_precision),' \sigma(3)=',num2str(sd_G5_3,str_precision)]};
            str_4={[str_Y,'(4)=',num2str(V_G5_4,str_precision),' r(4)=',num2str(r_G5_4,str_precision),' \sigma(4)=',num2str(sd_G5_4,str_precision)]};
            str_5={[str_Y,'(5)=',num2str(V_G5_5,str_precision),' r(5)=',num2str(r_G5_5,str_precision),' \sigma(5)=',num2str(sd_G5_5,str_precision)]};
            h_1=semilogx(X_interp,G5_1,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_2=semilogx(X_interp,G5_2,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_3=semilogx(X_interp,G5_3,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_4=semilogx(X_interp,G5_4,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_5=semilogx(X_interp,G5_5,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            semilogx(X_interp,y_G5_U,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on
            semilogx(X_interp,y_G5_L,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on           
            line([exp(B1) exp(B1)],[0 max(G5_1)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2) exp(B2)],[0 max(G5_2)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3) exp(B3)],[0 max(G5_3)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B4) exp(B4)],[0 max(G5_4)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B5) exp(B5)],[0 max(G5_5)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B1-sqrt(2*log(2))*sd_G5_1) exp(B1+sqrt(2*log(2))*sd_G5_1)],[max(G5_1)/2 max(G5_1)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2-sqrt(2*log(2))*sd_G5_2) exp(B2+sqrt(2*log(2))*sd_G5_2)],[max(G5_2)/2 max(G5_2)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3-sqrt(2*log(2))*sd_G5_3) exp(B3+sqrt(2*log(2))*sd_G5_3)],[max(G5_3)/2 max(G5_3)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B4-sqrt(2*log(2))*sd_G5_4) exp(B4+sqrt(2*log(2))*sd_G5_4)],[max(G5_4)/2 max(G5_4)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B5-sqrt(2*log(2))*sd_G5_5) exp(B5+sqrt(2*log(2))*sd_G5_5)],[max(G5_5)/2 max(G5_5)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey);   
            xlabel(str_x,'FontSize',font_size,'FontWeight','bold');ylabel(str_y,'FontSize',font_size,'FontWeight','bold');
            title(['5-modal: ','b=',num2str(b_5,str_precision),' s=',num2str(rmse_5,str_precision),' R^2=',num2str(R2adj_5,str_precision)],'FontSize',font_size,'FontWeight','bold');
            legendText=[];  
            legendText=[legendText,str_G5,str_1,str_2,str_3,str_4,str_5];
            legend([h_G5,h_1,h_2,h_3,h_4,h_5],legendText(:),'Location','EastOutside','FontSize',font_size);
            ylimits = get(gca,'ylim'); set(gca,'ylim',[0,ylimits(2)]);
            set(gcf,'PaperPositionMode','auto');  
            set(gca,'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XColor','k','YColor','k','LineWidth',line_width,'FontSize',font_size);
            if isequal(FLAG_plot_save,1)
                file_name = fullfile(path_save,['GMM5']); 
                print('-djpeg','-r200',file_name);                                                       
                close              
            end
                        
            if isequal(FLAG_plot_visible,0)                 
                figure; set(gcf, 'color','white', 'visible','off','units','normalized','outerposition',[0 0 1 1]);             
            else
                figure; set(gcf, 'color','white', 'visible','on','units','normalized','outerposition',[0 0 1 1]);                 
            end;
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
            hold on       
            h_G6 = semilogx(X_interp,y_G6,'color',bright,'LineWidth',2*line_width);
            hold on
            A1=fitresult6.a1;
            A2=fitresult6.a2;
            A3=fitresult6.a3;
            A4=fitresult6.a4;
            A5=fitresult6.a5;
            A6=fitresult6.a6;
            B1=fitresult6.b1;
            B2=fitresult6.b2;
            B3=fitresult6.b3;
            B4=fitresult6.b4;
            B5=fitresult6.b5;
            B6=fitresult6.b6;
            C1=fitresult6.c1;
            C2=fitresult6.c2;
            C3=fitresult6.c3;
            C4=fitresult6.c4;
            C5=fitresult6.c5;
            C6=fitresult6.c6;
            G6_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_5 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+A5*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            G6_6 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+A6*exp(-((logX-B6)/C6).^2); % 6-modal log-normal distribution
            V_G6_1 = trapz(logX,G6_1);
            V_G6_2 = trapz(logX,G6_2);
            V_G6_3 = trapz(logX,G6_3);
            V_G6_4 = trapz(logX,G6_4);
            V_G6_5 = trapz(logX,G6_5);
            V_G6_6 = trapz(logX,G6_6);
            r_G6_1 = exp(B1);
            r_G6_2 = exp(B2);
            r_G6_3 = exp(B3);
            r_G6_4 = exp(B4);
            r_G6_5 = exp(B5);
            r_G6_6 = exp(B6);
            sd_G6_1 = (1/sqrt(2))*(C1);    
            sd_G6_2 = (1/sqrt(2))*(C2);    
            sd_G6_3 = (1/sqrt(2))*(C3);    
            sd_G6_4 = (1/sqrt(2))*(C4);    
            sd_G6_5 = (1/sqrt(2))*(C5);    
            sd_G6_6 = (1/sqrt(2))*(C6);    
            x=Y_interp; y=y_G6;  n_coeff=18;
            [res_6,sse_6,R2_6,rmse_6]=my_stats(x,y); b_6=sum(x-y)/n_interp; df_6=n_interp-n_coeff; R2adj_6=1-(1-R2_6)*(n_interp-1)/(df_6-1);              
            str_1={[str_Y,'(1)=',num2str(V_G6_1,str_precision),' r(1)=',num2str(r_G6_1,str_precision),' \sigma(1)=',num2str(sd_G6_1,str_precision)]};
            str_2={[str_Y,'(2)=',num2str(V_G6_2,str_precision),' r(2)=',num2str(r_G6_2,str_precision),' \sigma(2)=',num2str(sd_G6_2,str_precision)]}; 
            str_3={[str_Y,'(3)=',num2str(V_G6_3,str_precision),' r(3)=',num2str(r_G6_3,str_precision),' \sigma(3)=',num2str(sd_G6_3,str_precision)]};
            str_4={[str_Y,'(4)=',num2str(V_G6_4,str_precision),' r(4)=',num2str(r_G6_4,str_precision),' \sigma(4)=',num2str(sd_G6_4,str_precision)]};
            str_5={[str_Y,'(5)=',num2str(V_G6_5,str_precision),' r(5)=',num2str(r_G6_5,str_precision),' \sigma(5)=',num2str(sd_G6_5,str_precision)]};
            str_6={[str_Y,'(6)=',num2str(V_G6_6,str_precision),' r(6)=',num2str(r_G6_6,str_precision),' \sigma(6)=',num2str(sd_G6_6,str_precision)]};
            h_1=semilogx(X_interp,G6_1,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_2=semilogx(X_interp,G6_2,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_3=semilogx(X_interp,G6_3,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_4=semilogx(X_interp,G6_4,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_5=semilogx(X_interp,G6_5,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_6=semilogx(X_interp,G6_6,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            semilogx(X_interp,y_G6_U,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on
            semilogx(X_interp,y_G6_L,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on           
            line([exp(B1) exp(B1)],[0 max(G6_1)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2) exp(B2)],[0 max(G6_2)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3) exp(B3)],[0 max(G6_3)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B4) exp(B4)],[0 max(G6_4)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B5) exp(B5)],[0 max(G6_5)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B6) exp(B6)],[0 max(G6_6)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B1-sqrt(2*log(2))*sd_G6_1) exp(B1+sqrt(2*log(2))*sd_G6_1)],[max(G6_1)/2 max(G6_1)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2-sqrt(2*log(2))*sd_G6_2) exp(B2+sqrt(2*log(2))*sd_G6_2)],[max(G6_2)/2 max(G6_2)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3-sqrt(2*log(2))*sd_G6_3) exp(B3+sqrt(2*log(2))*sd_G6_3)],[max(G6_3)/2 max(G6_3)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B4-sqrt(2*log(2))*sd_G6_4) exp(B4+sqrt(2*log(2))*sd_G6_4)],[max(G6_4)/2 max(G6_4)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B5-sqrt(2*log(2))*sd_G6_5) exp(B5+sqrt(2*log(2))*sd_G6_5)],[max(G6_5)/2 max(G6_5)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B6-sqrt(2*log(2))*sd_G6_6) exp(B6+sqrt(2*log(2))*sd_G6_6)],[max(G6_6)/2 max(G6_6)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey);   
            xlabel(str_x,'FontSize',font_size,'FontWeight','bold');ylabel(str_y,'FontSize',font_size,'FontWeight','bold');
            title(['6-modal: ','b=',num2str(b_6,str_precision),' s=',num2str(rmse_6,str_precision),' R^2=',num2str(R2adj_6,str_precision)],'FontSize',font_size,'FontWeight','bold');
            legendText=[];  
            legendText=[legendText,str_G6,str_1,str_2,str_3,str_4,str_5,str_6];
            legend([h_G6,h_1,h_2,h_3,h_4,h_5,h_6],legendText(:),'Location','EastOutside','FontSize',font_size);      
            ylimits = get(gca,'ylim'); set(gca,'ylim',[0,ylimits(2)]);
            set(gcf,'PaperPositionMode','auto');  
            set(gca,'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XColor','k','YColor','k','LineWidth',line_width,'FontSize',font_size);
            if isequal(FLAG_plot_save,1)
                file_name = fullfile(path_save,['GMM6']); 
                print('-djpeg','-r200',file_name);                                                       
                close              
            end

            if isequal(FLAG_plot_visible,0)
                figure; set(gcf, 'color','white', 'visible','off','units','normalized','outerposition',[0 0 1 1]);
            else   
                figure; set(gcf, 'color','white', 'visible','on','units','normalized','outerposition',[0 0 1 1]);    
            end        
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
            hold on      
            h_G7 = semilogx(X_interp,y_G7,'color',bright,'LineWidth',2*line_width);
            hold on
            A1=fitresult7.a1;
            A2=fitresult7.a2;
            A3=fitresult7.a3;
            A4=fitresult7.a4;
            A5=fitresult7.a5;
            A6=fitresult7.a6;
            A7=fitresult7.a7;
            B1=fitresult7.b1;
            B2=fitresult7.b2;
            B3=fitresult7.b3;
            B4=fitresult7.b4;
            B5=fitresult7.b5;
            B6=fitresult7.b6;
            B7=fitresult7.b7;
            C1=fitresult7.c1;
            C2=fitresult7.c2;
            C3=fitresult7.c3;
            C4=fitresult7.c4;
            C5=fitresult7.c5;
            C6=fitresult7.c6;
            C7=fitresult7.c7;
            G7_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_5 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+A5*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_6 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+A6*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            G7_7 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+A7*exp(-((logX-B7)/C7).^2); % 7-modal log-normal distribution
            V_G7_1 = trapz(logX,G7_1);
            V_G7_2 = trapz(logX,G7_2);
            V_G7_3 = trapz(logX,G7_3);
            V_G7_4 = trapz(logX,G7_4);
            V_G7_5 = trapz(logX,G7_5);
            V_G7_6 = trapz(logX,G7_6);
            V_G7_7 = trapz(logX,G7_7);
            r_G7_1 = exp(B1);
            r_G7_2 = exp(B2);
            r_G7_3 = exp(B3);
            r_G7_4 = exp(B4);
            r_G7_5 = exp(B5);
            r_G7_6 = exp(B6);
            r_G7_7 = exp(B7);
            sd_G7_1 = (1/sqrt(2))*(C1);    
            sd_G7_2 = (1/sqrt(2))*(C2);    
            sd_G7_3 = (1/sqrt(2))*(C3);    
            sd_G7_4 = (1/sqrt(2))*(C4);    
            sd_G7_5 = (1/sqrt(2))*(C5);    
            sd_G7_6 = (1/sqrt(2))*(C6);    
            sd_G7_7 = (1/sqrt(2))*(C7);    
            x=Y_interp; y=y_G7;  n_coeff=21;
            [res_7,sse_7,R2_7,rmse_7]=my_stats(x,y); b_7=sum(x-y)/n_interp; df_7=n_interp-n_coeff; R2adj_7=1-(1-R2_7)*(n_interp-1)/(df_7-1);              
            str_1={[str_Y,'(1)=',num2str(V_G7_1,str_precision),' r(1)=',num2str(r_G7_1,str_precision),' \sigma(1)=',num2str(sd_G7_1,str_precision)]};
            str_2={[str_Y,'(2)=',num2str(V_G7_2,str_precision),' r(2)=',num2str(r_G7_2,str_precision),' \sigma(2)=',num2str(sd_G7_2,str_precision)]}; 
            str_3={[str_Y,'(3)=',num2str(V_G7_3,str_precision),' r(3)=',num2str(r_G7_3,str_precision),' \sigma(3)=',num2str(sd_G7_3,str_precision)]};
            str_4={[str_Y,'(4)=',num2str(V_G7_4,str_precision),' r(4)=',num2str(r_G7_4,str_precision),' \sigma(4)=',num2str(sd_G7_4,str_precision)]};
            str_5={[str_Y,'(5)=',num2str(V_G7_5,str_precision),' r(5)=',num2str(r_G7_5,str_precision),' \sigma(5)=',num2str(sd_G7_5,str_precision)]};
            str_6={[str_Y,'(6)=',num2str(V_G7_6,str_precision),' r(6)=',num2str(r_G7_6,str_precision),' \sigma(6)=',num2str(sd_G7_6,str_precision)]};
            str_7={[str_Y,'(7)=',num2str(V_G7_7,str_precision),' r(7)=',num2str(r_G7_7,str_precision),' \sigma(7)=',num2str(sd_G7_7,str_precision)]};
            h_1=semilogx(X_interp,G7_1,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_2=semilogx(X_interp,G7_2,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_3=semilogx(X_interp,G7_3,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_4=semilogx(X_interp,G7_4,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_5=semilogx(X_interp,G7_5,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_6=semilogx(X_interp,G7_6,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            h_7=semilogx(X_interp,G7_7,'color',mid,'LineWidth',line_width,'LineStyle','-');
            hold on
            semilogx(X_interp,y_G7_U,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on
            semilogx(X_interp,y_G7_L,'color','k','LineWidth',line_width,'LineStyle','--');
            hold on           
            line([exp(B1) exp(B1)],[0 max(G7_1)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2) exp(B2)],[0 max(G7_2)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3) exp(B3)],[0 max(G7_3)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B4) exp(B4)],[0 max(G7_4)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B5) exp(B5)],[0 max(G7_5)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B6) exp(B6)],[0 max(G7_6)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B7) exp(B7)],[0 max(G7_7)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B1-sqrt(2*log(2))*sd_G7_1) exp(B1+sqrt(2*log(2))*sd_G7_1)],[max(G7_1)/2 max(G7_1)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B2-sqrt(2*log(2))*sd_G7_2) exp(B2+sqrt(2*log(2))*sd_G7_2)],[max(G7_2)/2 max(G7_2)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B3-sqrt(2*log(2))*sd_G7_3) exp(B3+sqrt(2*log(2))*sd_G7_3)],[max(G7_3)/2 max(G7_3)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B4-sqrt(2*log(2))*sd_G7_4) exp(B4+sqrt(2*log(2))*sd_G7_4)],[max(G7_4)/2 max(G7_4)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B5-sqrt(2*log(2))*sd_G7_5) exp(B5+sqrt(2*log(2))*sd_G7_5)],[max(G7_5)/2 max(G7_5)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B6-sqrt(2*log(2))*sd_G7_6) exp(B6+sqrt(2*log(2))*sd_G7_6)],[max(G7_6)/2 max(G7_6)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            line([exp(B7-sqrt(2*log(2))*sd_G7_7) exp(B7+sqrt(2*log(2))*sd_G7_7)],[max(G7_7)/2 max(G7_7)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
            hold on
            h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey);   
            xlabel(str_x,'FontSize',font_size,'FontWeight','bold');ylabel(str_y,'FontSize',font_size,'FontWeight','bold');
            title(['7-modal: ','b=',num2str(b_7,str_precision),' s=',num2str(rmse_7,str_precision),' R^2=',num2str(R2adj_7,str_precision)],'FontSize',font_size,'FontWeight','bold');
            legendText=[];  
            legendText=[legendText,str_G7,str_1,str_2,str_3,str_4,str_5,str_6,str_7];
            legend([h_G7,h_1,h_2,h_3,h_4,h_5,h_6,h_7],legendText(:),'Location','EastOutside','FontSize',font_size);
            ylimits = get(gca,'ylim'); set(gca,'ylim',[0,ylimits(2)]);
            set(gcf,'PaperPositionMode','auto');  
            set(gca,'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XColor','k','YColor','k','LineWidth',line_width,'FontSize',font_size);
            if isequal(FLAG_plot_save,1)
                file_name = fullfile(path_save,['GMM7']); 
                print('-djpeg','-r200',file_name);                                                       
                close              
            end
            
            if isequal(FLAG_plot_visible,0)                 
                figure; set(gcf, 'color','white', 'visible','off','units','normalized','outerposition',[0 0 1 1]);             
            else
                figure; set(gcf, 'color','white', 'visible','on','units','normalized','outerposition',[0 0 1 1]);                 
            end;
            if N_interp>24
                h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey); 
                hold on      
                h_G8 = semilogx(X_interp,y_G8,'color',bright,'LineWidth',2*line_width);
                hold on     
                A1=fitresult8.a1;
                A2=fitresult8.a2;
                A3=fitresult8.a3;
                A4=fitresult8.a4;
                A5=fitresult8.a5;
                A6=fitresult8.a6;
                A7=fitresult8.a7;
                A8=fitresult8.a8;
                B1=fitresult8.b1;
                B2=fitresult8.b2;
                B3=fitresult8.b3;
                B4=fitresult8.b4;
                B5=fitresult8.b5;
                B6=fitresult8.b6;
                B7=fitresult8.b7;
                B8=fitresult8.b8;
                C1=fitresult8.c1;
                C2=fitresult8.c2;
                C3=fitresult8.c3;
                C4=fitresult8.c4;
                C5=fitresult8.c5;
                C6=fitresult8.c6;
                C7=fitresult8.c7;
                C8=fitresult8.c8;
                G8_1 = A1*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
                G8_2 = 0*exp(-((logX-B1)/C1).^2)+A2*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
                G8_3 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+A3*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
                G8_4 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+A4*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
                G8_5 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+A5*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
                G8_6 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+A6*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
                G8_7 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+A7*exp(-((logX-B7)/C7).^2)+0*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
                G8_8 = 0*exp(-((logX-B1)/C1).^2)+0*exp(-((logX-B2)/C2).^2)+0*exp(-((logX-B3)/C3).^2)+0*exp(-((logX-B4)/C4).^2)+0*exp(-((logX-B5)/C5).^2)+0*exp(-((logX-B6)/C6).^2)+0*exp(-((logX-B7)/C7).^2)+A8*exp(-((logX-B8)/C8).^2); % 8-modal log-normal distribution
                V_G8_1 = trapz(logX,G8_1);
                V_G8_2 = trapz(logX,G8_2);
                V_G8_3 = trapz(logX,G8_3);
                V_G8_4 = trapz(logX,G8_4);
                V_G8_5 = trapz(logX,G8_5);
                V_G8_6 = trapz(logX,G8_6);
                V_G8_7 = trapz(logX,G8_7);
                V_G8_8 = trapz(logX,G8_8);
                r_G8_1 = exp(B1);
                r_G8_2 = exp(B2);
                r_G8_3 = exp(B3);
                r_G8_4 = exp(B4);
                r_G8_5 = exp(B5);
                r_G8_6 = exp(B6);
                r_G8_7 = exp(B7);
                r_G8_8 = exp(B8);
                sd_G8_1 = (1/sqrt(2))*(C1);    
                sd_G8_2 = (1/sqrt(2))*(C2);    
                sd_G8_3 = (1/sqrt(2))*(C3);    
                sd_G8_4 = (1/sqrt(2))*(C4);    
                sd_G8_5 = (1/sqrt(2))*(C5);    
                sd_G8_6 = (1/sqrt(2))*(C6);    
                sd_G8_7 = (1/sqrt(2))*(C7);    
                sd_G8_8 = (1/sqrt(2))*(C8);    
                x=Y_interp; y=y_G8;  n_coeff=24;
                [res_8,sse_8,R2_8,rmse_8]=my_stats(x,y); b_8=sum(x-y)/n_interp; df_8=n_interp-n_coeff; R2adj_8=1-(1-R2_8)*(n_interp-1)/(df_8-1);              
                str_1={[str_Y,'(1)=',num2str(V_G8_1,str_precision),' r(1)=',num2str(r_G8_1,str_precision),' \sigma(1)=',num2str(sd_G8_1,str_precision)]};
                str_2={[str_Y,'(2)=',num2str(V_G8_2,str_precision),' r(2)=',num2str(r_G8_2,str_precision),' \sigma(2)=',num2str(sd_G8_2,str_precision)]}; 
                str_3={[str_Y,'(3)=',num2str(V_G8_3,str_precision),' r(3)=',num2str(r_G8_3,str_precision),' \sigma(3)=',num2str(sd_G8_3,str_precision)]};
                str_4={[str_Y,'(4)=',num2str(V_G8_4,str_precision),' r(4)=',num2str(r_G8_4,str_precision),' \sigma(4)=',num2str(sd_G8_4,str_precision)]};
                str_5={[str_Y,'(5)=',num2str(V_G8_5,str_precision),' r(5)=',num2str(r_G8_5,str_precision),' \sigma(5)=',num2str(sd_G8_5,str_precision)]};
                str_6={[str_Y,'(6)=',num2str(V_G8_6,str_precision),' r(6)=',num2str(r_G8_6,str_precision),' \sigma(6)=',num2str(sd_G8_6,str_precision)]};
                str_7={[str_Y,'(7)=',num2str(V_G8_7,str_precision),' r(7)=',num2str(r_G8_7,str_precision),' \sigma(7)=',num2str(sd_G8_7,str_precision)]};
                str_8={[str_Y,'(8)=',num2str(V_G8_8,str_precision),' r(8)=',num2str(r_G8_8,str_precision),' \sigma(8)=',num2str(sd_G8_8,str_precision)]};
                h_1=semilogx(X_interp,G8_1,'color',mid,'LineWidth',line_width,'LineStyle','-');
                hold on
                h_2=semilogx(X_interp,G8_2,'color',mid,'LineWidth',line_width,'LineStyle','-');
                hold on
                h_3=semilogx(X_interp,G8_3,'color',mid,'LineWidth',line_width,'LineStyle','-');
                hold on
                h_4=semilogx(X_interp,G8_4,'color',mid,'LineWidth',line_width,'LineStyle','-');
                hold on
                h_5=semilogx(X_interp,G8_5,'color',mid,'LineWidth',line_width,'LineStyle','-');
                hold on
                h_6=semilogx(X_interp,G8_6,'color',mid,'LineWidth',line_width,'LineStyle','-');
                hold on
                h_7=semilogx(X_interp,G8_7,'color',mid,'LineWidth',line_width,'LineStyle','-');
                hold on
                h_8=semilogx(X_interp,G8_8,'color',mid,'LineWidth',line_width,'LineStyle','-');
                hold on
                semilogx(X_interp,y_G8_U,'color','k','LineWidth',line_width,'LineStyle','--');
                hold on
                semilogx(X_interp,y_G8_L,'color','k','LineWidth',line_width,'LineStyle','--');
                hold on           
                line([exp(B1) exp(B1)],[0 max(G8_1)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B2) exp(B2)],[0 max(G8_2)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B3) exp(B3)],[0 max(G8_3)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B4) exp(B4)],[0 max(G8_4)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B5) exp(B5)],[0 max(G8_5)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B6) exp(B6)],[0 max(G8_6)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B7) exp(B7)],[0 max(G8_7)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B8) exp(B8)],[0 max(G8_8)],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B1-sqrt(2*log(2))*sd_G8_1) exp(B1+sqrt(2*log(2))*sd_G8_1)],[max(G8_1)/2 max(G8_1)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B2-sqrt(2*log(2))*sd_G8_2) exp(B2+sqrt(2*log(2))*sd_G8_2)],[max(G8_2)/2 max(G8_2)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B3-sqrt(2*log(2))*sd_G8_3) exp(B3+sqrt(2*log(2))*sd_G8_3)],[max(G8_3)/2 max(G8_3)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B4-sqrt(2*log(2))*sd_G8_4) exp(B4+sqrt(2*log(2))*sd_G8_4)],[max(G8_4)/2 max(G8_4)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B5-sqrt(2*log(2))*sd_G8_5) exp(B5+sqrt(2*log(2))*sd_G8_5)],[max(G8_5)/2 max(G8_5)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B6-sqrt(2*log(2))*sd_G8_6) exp(B6+sqrt(2*log(2))*sd_G8_6)],[max(G8_6)/2 max(G8_6)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B7-sqrt(2*log(2))*sd_G8_7) exp(B7+sqrt(2*log(2))*sd_G8_7)],[max(G8_7)/2 max(G8_7)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                line([exp(B8-sqrt(2*log(2))*sd_G8_8) exp(B8+sqrt(2*log(2))*sd_G8_8)],[max(G8_8)/2 max(G8_8)/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');
                hold on
                h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey);   
                xlabel(str_x,'FontSize',font_size,'FontWeight','bold');ylabel(str_y,'FontSize',font_size,'FontWeight','bold');
                title(['8-modal: ','b=',num2str(b_8,str_precision),' s=',num2str(rmse_8,str_precision),' R^2=',num2str(R2adj_8,str_precision)],'FontSize',font_size,'FontWeight','bold');
                legendText=[];  
                legendText=[legendText,str_G8,str_1,str_2,str_3,str_4,str_5,str_6,str_7,str_8];
                legend([h_G8,h_1,h_2,h_3,h_4,h_5,h_6,h_7,h_8],legendText(:),'Location','EastOutside','FontSize',font_size);      
            else
                set(gcf,'Visible','off');
            end
            ylimits = get(gca,'ylim'); set(gca,'ylim',[0,ylimits(2)]);
            set(gcf,'PaperPositionMode','auto');  
            set(gca,'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XColor','k','YColor','k','LineWidth',line_width,'FontSize',font_size);
            if isequal(FLAG_plot_save,1)
                file_name = fullfile(path_save,['GMM8']); 
                print('-djpeg','-r200',file_name);                                                       
                close              
            end   
        end

    %% STORE FIT PARAMETERS
    GMM_proc=[];    
    GMM_proc=[GMM_proc;zeros(1,30)]; 
    GMM_proc=[GMM_proc;zeros(1,30)]; 
    GMM_proc=[GMM_proc;zeros(1,30)]; 
    GMM_proc=[GMM_proc;zeros(1,30)];     
    temp1=[V_G1_1,ones(1,7)*NaN;r_G1_1,ones(1,7)*NaN;sd_G1_1,ones(1,7)*NaN]; 
    tempA=sortrows(temp1',2); tempB=[tempA(:,1)',tempA(:,2)',tempA(:,3)'];    
    GMM_proc=[GMM_proc;NaN,1,double(b_1),rmse_1,double(R2adj_1),V_G1,tempB]; %1 mode    
    temp2=[V_G2_1,V_G2_2,ones(1,6)*NaN;r_G2_1,r_G2_2,ones(1,6)*NaN;sd_G2_1,sd_G2_2,ones(1,6)*NaN]; 
    tempA=sortrows(temp2',2); tempB=[tempA(:,1)',tempA(:,2)',tempA(:,3)'];
    GMM_proc=[GMM_proc;NaN,2,double(b_2),rmse_2,double(R2adj_2),V_G2,tempB]; %2 modes    
    temp3=[V_G3_1,V_G3_2,V_G3_3,ones(1,5)*NaN;r_G3_1,r_G3_2,r_G3_3,ones(1,5)*NaN;sd_G3_1,sd_G3_2,sd_G3_3,ones(1,5)*NaN];                                
    tempA=sortrows(temp3',2); tempB=[tempA(:,1)',tempA(:,2)',tempA(:,3)'];    
    GMM_proc=[GMM_proc;NaN,3,double(b_3),rmse_3,double(R2adj_3),V_G3,tempB]; %3 modes    
    temp4=[V_G4_1,V_G4_2,V_G4_3,V_G4_4,ones(1,4)*NaN;r_G4_1,r_G4_2,r_G4_3,r_G4_4,ones(1,4)*NaN;sd_G4_1,sd_G4_2,sd_G4_3,sd_G4_4,ones(1,4)*NaN];                                
    tempA=sortrows(temp4',2); tempB=[tempA(:,1)',tempA(:,2)',tempA(:,3)'];    
    GMM_proc=[GMM_proc;NaN,4,double(b_4),rmse_4,double(R2adj_4),V_G4,tempB]; %4 modes    
    temp5=[V_G5_1,V_G5_2,V_G5_3,V_G5_4,V_G5_5,ones(1,3)*NaN;r_G5_1,r_G5_2,r_G5_3,r_G5_4,r_G5_5,ones(1,3)*NaN;sd_G5_1,sd_G5_2,sd_G5_3,sd_G5_4,sd_G5_5,ones(1,3)*NaN];                                
    tempA=sortrows(temp5',2); tempB=[tempA(:,1)',tempA(:,2)',tempA(:,3)'];    
    GMM_proc=[GMM_proc;NaN,5,double(b_5),rmse_5,double(R2adj_5),V_G5,tempB]; %5 modes    
    temp6=[V_G6_1,V_G6_2,V_G6_3,V_G6_4,V_G6_5,V_G6_6,ones(1,2)*NaN;r_G6_1,r_G6_2,r_G6_3,r_G6_4,r_G6_5,r_G6_6,ones(1,2)*NaN;sd_G6_1,sd_G6_2,sd_G6_3,sd_G6_4,sd_G6_5,sd_G6_6,ones(1,2)*NaN];                                
    tempA=sortrows(temp6',2); tempB=[tempA(:,1)',tempA(:,2)',tempA(:,3)'];    
    GMM_proc=[GMM_proc;NaN,6,double(b_6),rmse_6,double(R2adj_6),V_G6,tempB]; %6 modes    
    temp7=[V_G7_1,V_G7_2,V_G7_3,V_G7_4,V_G7_5,V_G7_6,V_G7_7,ones(1,1)*NaN;r_G7_1,r_G7_2,r_G7_3,r_G7_4,r_G7_5,r_G7_6,r_G7_7,ones(1,1)*NaN;sd_G7_1,sd_G7_2,sd_G7_3,sd_G7_4,sd_G7_5,sd_G7_6,sd_G7_7,ones(1,1)*NaN];                                
    tempA=sortrows(temp7',2); tempB=[tempA(:,1)',tempA(:,2)',tempA(:,3)'];    
    GMM_proc=[GMM_proc;NaN,7,double(b_7),rmse_7,double(R2adj_7),V_G7,tempB]; %7 modes    
    if N_interp>24
        temp8=[V_G8_1,V_G8_2,V_G8_3,V_G8_4,V_G8_5,V_G8_6,V_G8_7,V_G8_8;r_G8_1,r_G8_2,r_G8_3,r_G8_4,r_G8_5,r_G8_6,r_G8_7,r_G8_8;sd_G8_1,sd_G8_2,sd_G8_3,sd_G8_4,sd_G8_5,sd_G8_6,sd_G8_7,sd_G8_8];                                                                           
        tempA=sortrows(temp8',2); tempB=[tempA(:,1)',tempA(:,2)',tempA(:,3)'];    
        GMM_proc=[GMM_proc;NaN,8,double(b_8),rmse_8,double(R2adj_8),V_G8,tempB]; %8 modes    
    else
        GMM_proc=[GMM_proc;NaN,8,ones(1,28)*NaN]; %8 modes 
    end
    GMM_proc(isnan(GMM_proc))=0;       

    GMM_test=zeros(8,18);       
    if isequal(FLAG_autodetect,0)  
        prompt={'Enter the optimal number of modes:'};
        name='Manual n-modes selection';
        numlines=1;
        defaultanswer={'3'};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        n_modes=str2double(answer);
    else
        % BEST GMM MODEL FIT                
        V_diff=abs([GMM_proc(5:12,6)]-V_trapz);
        opt_V=find(V_diff==min(V_diff));    
        opt_s=find(GMM_proc(:,4)==min(GMM_proc(5:12,4)))-1;    
        opt_R2=find(GMM_proc(:,5)==max(GMM_proc(5:12,5)))-1;        
        F0=1.96;
        N=n_interp;
        R=[];
        FLAG_stop=0;
        for j=1:8       
            if isequal(FLAG_autodetect,1)
                if isequal(FLAG_method,0)                
                    % MULTIPLE REGRESSION METHOD
                    if j>1    
                        R2adj_F=GMM_proc(j+4,5); % n
                        R2adj_R=GMM_proc(j+4-1,5); % n-1                  
                        F=((R2adj_F^2-R2adj_R)/3) / ((1-R2adj_F)/(N-(3*j)-1));
                        % alpha F                   alpha       F
                        % 0.80	1.281551565545		0.999       3.290526731492
                        % 0.90	1.644853626951		0.9999      3.890591886413
                        % 0.95	1.959963984540		0.99999     4.417173413469
                        % 0.98	2.326347874041		0.999999	4.891638475699
                        % 0.99	2.575829303549		0.9999999	5.326723886384
                        % 0.995	2.807033768344		0.99999999	5.730728868236
                        % 0.998	3.090232306168		0.999999999	6.109410204869                
                        if F<F0
                            significant=1;
                            n_modes=j;
                        else                        
                            significant=0;                                                           
                            disp(['n=',num2str(j),': Non-significant change in R2']);                               
                        end
                        GMM_test(j,:)=[j,GMM_proc(j+4-1,5),GMM_proc(j+4,5),ones(1,10),F,F0,NaN,significant,j];
                    elseif isequal(j,1)
                        GMM_test(j,:)=[j,GMM_proc(j+4,5),NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1];
                    end
                end                
                if isequal(FLAG_method,1)
                    % FISHER TRANSFORM METHOD
                    if j>1                    
                        z1=0.5*log( (1+sqrt(GMM_proc(j+4-1,5))) / (1-sqrt(GMM_proc(j+4-1,5))) );
                        z2=0.5*log( (1+sqrt(GMM_proc(j+4,5))) / (1-sqrt(GMM_proc(j+4,5))) );
                        SE1=1/sqrt(N-3);
                        SE2=1/sqrt(N-3);
                        t_Welch=abs( (z1-z2)/sqrt( (1/(N-3))+(1/(N-3)) ));
                        % alpha F                   alpha       F
                        % 0.80	1.281551565545		0.999       3.290526731492
                        % 0.90	1.644853626951		0.9999      3.890591886413
                        % 0.95	1.959963984540		0.99999     4.417173413469
                        % 0.98	2.326347874041		0.999999	4.891638475699
                        % 0.99	2.575829303549		0.9999999	5.326723886384
                        % 0.995	2.807033768344		0.99999999	5.730728868236
                        % 0.998	3.090232306168		0.999999999	6.109410204869                
                        CI1_U=z1+(F0*SE1);
                        CI1_L=z1-(F0*SE1);
                        CI2_U=z2+(F0*SE2);
                        CI2_L=z2-(F0*SE2);                                               
                        if GMM_proc(j+4,5)>GMM_proc(j+4-1,5) % R2(n)>R2(n-1)
                            ordered=1;
                            if CI2_L>CI1_U
                                separated=1;
                                t_test=NaN;                                
                                significant=1;
                                n_modes=j;
                            else
                                separated=0;
                                if t_Welch>F0
                                    t_test=1;
                                    significant=1;    
                                    n_modes=j;
                                else
                                    t_test=0;
                                    significant=0;                                                                                                
                                    disp(['n=',num2str(j),': Non-significant change in R2']);                               
                                end       
                            end                            
                        else
                            ordered=0;
                            separated=NaN;                            
                            t_test=NaN;
                            significant=0;                                                                                        
                            disp(['n=',num2str(j),': Non-significant change in R2']);                              
                        end
                        GMM_test(j,:)=[j,GMM_proc(j+4-1,5),GMM_proc(j+4,5),ordered,z1,z2,SE1,SE2,CI1_L,CI1_U,CI2_L,CI2_U,separated,t_Welch,F0,t_test,significant,j];
                    elseif isequal(j,1)
                        GMM_test(j,:)=[j,GMM_proc(j+4,5),NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1];
                    end                
                end            
                if isequal(FLAG_method,2)                
                    % SUBJECTIVE METHOD
                    GMM_proc(j+2,3)=b(j)
                    GMM_proc(j+2,4)=s(j)
                    GMM_proc(j+2,5)=R2(j)
                    GMM_proc(j+2,6)=V_GMM(j)
                    GMM_proc(j+2,7:14)=V_GMM(1:8)                   
                    if j>1
                        if abs(GMM_proc(j+4-1,3))<abs(GMM_proc(j+4,3)) || GMM_proc(j+4-1,4)<GMM_proc(j+4,4) || GMM_proc(j+4-1,5)>GMM_proc(j+4,5)
                            significant=0;                                                           
                            disp(['n=',num2str(j),': Non-significant change in R2']);      
                        else            
                           significant=1;
                           n_modes=j;
                        end  
                        GMM_test(j,:)=[j,GMM_proc(j+4-1,5),GMM_proc(j+4,5),ones(1,10),NaN,NaN,NaN,significant,j];                        
                    elseif isequal(j,1)
                        GMM_test(j,:)=[j,GMM_proc(j+4,5),NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1];
                    end                       
                end
            else
            end
        end        
        stop_row1=find(GMM_test(:,4)==0);
        stop_row2=find(GMM_test(:,17)==0);
        if isequal(isempty(stop_row1),1)
            if isequal(isempty(stop_row2),1)
                n_modes=8;
            else
                n_modes=stop_row2(1)-1;
            end
        else
            if isequal(isempty(stop_row2),1)
                n_modes=stop_row1(1)-1;
            else
                if stop_row1(1)<stop_row2(1)
                    n_modes=stop_row1(1)-1;
                else                
                    n_modes=stop_row2(1)-1;
                end   
            end
        end    
        for k=2:8       
            if k>n_modes        
                GMM_test(k,18)=n_modes;            
            end
        end    
    end   
    if isequal(n_modes,1)
        y_best=y_G1;
        y_best_U=y_G1_U;
        y_best_L=y_G1_L;
        V_GMM=V_G1; 
        b_GMM=b_1;
        s_GMM=rmse_1;
        R2_GMM=R2adj_1;
    elseif isequal(n_modes,2)         
        y_best=y_G2;         
        y_best_U=y_G2_U;
        y_best_L=y_G2_L;
        V_GMM=V_G2;   
        b_GMM=b_2;
        s_GMM=rmse_2;
        R2_GMM=R2adj_2;
    elseif isequal(n_modes,3)         
        y_best=y_G3;         
        y_best_U=y_G3_U;
        y_best_L=y_G3_L;
        V_GMM=V_G3;                    
        b_GMM=b_3;
        s_GMM=rmse_3;
        R2_GMM=R2adj_3;
    elseif isequal(n_modes,4)         
        y_best=y_G4;         
        y_best_U=y_G4_U;
        y_best_L=y_G4_L;
        V_GMM=V_G4;                    
        b_GMM=b_4;
        s_GMM=rmse_4;
        R2_GMM=R2adj_4;
    elseif isequal(n_modes,5)         
        y_best=y_G5;         
        y_best_U=y_G5_U;
        y_best_L=y_G5_L;
        V_GMM=V_G5;                    
        b_GMM=b_5;
        s_GMM=rmse_5;
        R2_GMM=R2adj_5;
    elseif isequal(n_modes,6)         
        y_best=y_G6;         
        V_GMM=V_G6;                    
        y_best_U=y_G6_U;
        y_best_L=y_G6_L;
        b_GMM=b_6;
        s_GMM=rmse_6;
        R2_GMM=R2adj_6;
    elseif isequal(n_modes,7)         
        y_best=y_G7;         
        y_best_U=y_G7_U;
        y_best_L=y_G7_L;
        V_GMM=V_G7;                    
        b_GMM=b_7;
        s_GMM=rmse_7;
        R2_GMM=R2adj_7;
    elseif isequal(n_modes,8)         
        y_best=y_G8;         
        y_best_U=y_G8_U;
        y_best_L=y_G8_L;
        V_GMM=V_G8;                    
        b_GMM=b_8;
        s_GMM=rmse_8;
        R2_GMM=R2adj_8;
    end        

    %% PLOT GMM BEST FIT versus AERONET (+ spline fit)
    if isequal(FLAG_plot_visible,0)        
        figure; set(gcf, 'color','white', 'visible','off','units','normalized','outerposition',[0 0 1 1]);        
    else        
        figure; set(gcf, 'color','white', 'visible','on','units','normalized','outerposition',[0 0 1 1]);        
    end    
    h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey);    
    hold on        
    legendText=[];                   
    x=Y_interp;  
    y=y_best;
    h_FIT=semilogx(X_interp,y,'Color',bright,'LineWidth',2*line_width); % LAR bi-modal fit    
    hold on    
    str_modes={};
    h_MODES={};
    for k=1:n_modes       
        y_mode=GMM_proc(4+n_modes,7+k-1)/(sqrt(2*pi)*GMM_proc(4+n_modes,23+k-1))*exp(-((logX-log(GMM_proc(4+n_modes,15+k-1)))/(sqrt(2)*GMM_proc(4+n_modes,23+k-1))).^2);
        h_MODES{k}=semilogx(X_interp,y_mode,'Color',mid,'LineWidth',line_width);
        hold on
        line([exp(log(GMM_proc(4+n_modes,15+k-1) )) exp(log(GMM_proc(4+n_modes,15+k-1)))],[0 GMM_proc(4+n_modes,7+k-1)/(sqrt(2*pi)*GMM_proc(4+n_modes,23+k-1))],'LineWidth',line_width,'Color',mid,'LineStyle','-');     
        hold on    
        line([exp(log(GMM_proc(4+n_modes,15+k-1) )-sqrt(2*log(2))*GMM_proc(4+n_modes,23+k-1) ) exp(log(GMM_proc(4+n_modes,15+k-1))+sqrt(2*log(2))*GMM_proc(4+n_modes,23+k-1))],[GMM_proc(4+n_modes,7+k-1)/(sqrt(2*pi)*GMM_proc(4+n_modes,23+k-1))/2 GMM_proc(4+n_modes,7+k-1)/(sqrt(2*pi)*GMM_proc(4+n_modes,23+k-1))/2],'LineWidth',line_width,'Color',mid,'LineStyle','-');     
        hold on   
        temp=strcat([str_Y,'(',num2str(k),')=',num2str(GMM_proc(4+n_modes,7+k-1),str_precision),' r(',num2str(k),')=',num2str(GMM_proc(4+n_modes,15+k-1),str_precision),' \sigma(',num2str(k),')=',num2str(GMM_proc(4+n_modes,23+k-1),str_precision)]); 
        str_modes = [str_modes,temp];
    end     
    h_INTERP=semilogx(X_interp,Y_interp,'Color','k','LineWidth',line_width);    
    hold on    
    h_AERONET=semilogx(X,Y,'s','MarkerEdgeColor','k','MarkerSize',marker_size,'MarkerFaceColor',grey);     
    xlabel(str_x,'FontSize',font_size,'FontWeight',font_weight);ylabel(str_y,'FontSize',font_size,'FontWeight',font_weight);    
    title(['GMM best fit: ','b=',num2str(b_GMM,str_precision),' s=',num2str(s_GMM,str_precision),' R^2=',num2str(R2_GMM,str_precision)],'FontSize',font_size,'FontWeight','bold');
    str_AERONET={['Raw data: ',str_Y,'=',num2str(V_trapz,str_precision),' [22 points]']}; 
    str_INTERP={[ 'Spline fit: ',str_Y,'=',num2str(V_interp,str_precision),' [',num2str(n_interp),' points]']}; 
    str_FIT = {[num2str(n_modes),'-modal fit:',str_Y,'=',num2str(V_GMM,str_precision)]};    
    legendText=[legendText,str_AERONET,str_INTERP,str_FIT,str_modes];    
    legend([h_AERONET;h_INTERP;h_FIT;h_MODES{1};h_MODES{2}],legendText(:),'Location','EastOutside','FontSize',font_size);      
    ylimits = get(gca,'ylim'); set(gca,'ylim',[0,ylimits(2)]);    
    set(gcf,'PaperPositionMode','auto');      
    set(gca,'Box','off','TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XColor','k','YColor','k','LineWidth',line_width,'FontSize',font_size);    
    if isequal(FLAG_plot_save,1)    
        file_name = fullfile(path_save,['BEST']);         
        print('-djpeg','-r200',file_name);                                                               
        close                      
    end   

    GMM_stats = [GMM_proc(5:12,2:6)];
    GMM_modes = [(1:8)',GMM_proc(5:12,7:30)];

    field1 = 'stats';  value1 = {GMM_stats};
    field2 = 'modes';  value2 = {GMM_modes};
    field3 = 'coef';  value3 = {GMM_coef};
    field4 = 'se';  value4 = {GMM_se};
    field5 = 'test';  value5 = {GMM_test};

    GMM = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);

    file_name = fullfile(strcat(path_save,strcat('RUN','.mat')));   
    save([file_name],'GMM_stats','GMM_modes','GMM_coef','GMM_se','GMM_test');    

    end
    
    function [res,sse,R2,rmse]=my_stats(x,y)

      % x=data;
      % y=estimate;

      [Q,R] = qr(x,0); % orthogonal triangular decomposition
      beta = R\(Q'*y);
      yhat = x*beta;
      res = y - yhat;  
      nobs = length(y);
      p = length(beta);
      dfe = nobs-p;
      dft = nobs-1;
      ybar = mean(y);
      sse = norm(res)^2;                    % sum of squared errors
      ssr = norm(yhat - ybar)^2;            % regression sum of squares
      sst = norm(y - ybar)^2;               % total sum of squares;
      mse = sse./dfe;                       % mean squared error
      rmse = sqrt(mse);                     % root mean squared error
      R2 = 1 - sse ./ sst;                  % R^2 value
      R2_adj = 1 - (sse./sst)*(dft./dfe);   % adjusted R^2

    end
    
                