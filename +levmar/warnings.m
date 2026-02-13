function warnings(flag)
ratio_reached_imax = mean(flag(:)==2);
if ratio_reached_imax
    warning([num2str(ratio_reached_imax*100,2) ' % reached iteration limit'])
end
ratio_not_converged = mean(flag(:)==1);
if ratio_not_converged
    warning([num2str(ratio_not_converged*100,2) ' % could not be reduced to error target'])
end