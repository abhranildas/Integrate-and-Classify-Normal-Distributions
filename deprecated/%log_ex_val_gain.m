function l=log_ex_val_gain(x,mu,v,prior,val_gain)
dim=length(mu);
maha=dot((x-mu'),(x-mu')/v,2);
halflogdetv=sum(log(diag(chol(v))));
l=-(maha+log(2*pi)*dim)/2-halflogdetv+log(prior*val_gain);
