function l=normal_lpr(x,mu_top,v_top,mu_bot,v_bot,prior_top)
maha_top=dot((x-mu_top),(x-mu_top)/v_top,2);
maha_bot=dot((x-mu_bot),(x-mu_bot)/v_bot,2);
[R_top,~]=chol(v_top);
[R_bot,~]=chol(v_bot);
logdet_top=2*sum(log(diag(R_top)));
logdet_bot=2*sum(log(diag(R_bot)));
l=(maha_bot-maha_top+logdet_bot-logdet_top)/2+log(prior_top/(1-prior_top));
