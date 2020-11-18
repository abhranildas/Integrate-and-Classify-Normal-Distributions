function l=bayes_dec_var(x,mu_top,v_top,mu_bot,v_bot,prior_top,vals)
maha_top=dot((x-mu_top),(x-mu_top)/v_top,2);
maha_bot=dot((x-mu_bot),(x-mu_bot)/v_bot,2);
logdet_top=2*sum(log(diag(chol(v_top))));
logdet_bot=2*sum(log(diag(chol(v_bot))));
l=(maha_bot-maha_top+logdet_bot-logdet_top)/2 ...
+log(prior_top/(1-prior_top))...
+log((vals(1,1)-vals(1,2))/(vals(2,2)-vals(2,1)));
