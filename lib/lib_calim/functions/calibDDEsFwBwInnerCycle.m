function [DDE_U]=calibDDEsFwBwUpdates(DDE_U,BCst,GradOp,param)

dde_theta_maxo =    param.dde_theta_maxo;
die_theta_mino =  1-param.die_theta_maxo;
die_theta_maxo =  1+param.die_theta_maxo;
if param.flag_time_smoothness
    dieTempralFreq  = param.dieTempralFreq;
else
    dieTempralFreq = ':';
end

gamma = 1.98/((param.nuo + 1.00001*param.Lips) );% step verified paper Chouzenoux et al. 2016

% forward-backward step
for itr = 1:param.JUo
    %% FW step: gradient step
    DDE_U  = DDE_U - gamma*(GradOp(DDE_U) - BCst) ;
    %% BW step: projection
    centralComponentVal = DDE_U(dieTempralFreq,param.dieSpatialFreq);
    gurc =  min(max(real(centralComponentVal),   die_theta_mino), die_theta_maxo);
    guic =  min(max(imag(centralComponentVal),  -dde_theta_maxo), dde_theta_maxo) ;
    zero_freq_component = gurc + 1i*guic;
    
    gur =  min(max(real(DDE_U), -dde_theta_maxo), dde_theta_maxo);
    gui =  min(max(imag(DDE_U), -dde_theta_maxo), dde_theta_maxo) ;
    DDE_U  =  gur + 1i*gui;
    
    DDE_U(dieTempralFreq,param.dieSpatialFreq) = zero_freq_component;
end
end

