function checkresiduals(rho_obs, rhodot_obs, rho_comp, rhodot_comp)
    res_rho = rho_obs-rho_comp;
    res_rhodot = rhodot_obs-rhodot_comp;
    
    average_rho_res = mean(res_rho)
    average_rhodot_res = mean(res_rhodot)

end
